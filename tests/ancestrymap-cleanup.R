library(LEA)
local({
    base <- if (dir.exists("ancestrymap/fixtures")) "ancestrymap" else "tests/ancestrymap"
    child <- normalizePath(file.path(base,"isolated-cleanup.R"))
    root <- tempfile("ancestrymap-isolated-"); dir.create(root)
    on.exit(unlink(root,recursive=TRUE))
    checks <- 0L
    execute <- function(format, text="a i 1\n", pattern=NULL, expected=NULL,
                        override=list(), mode="wrapper", missing=FALSE, output_dir=FALSE) {
        checks <<- checks + 1L
        work <- file.path(root,paste0("case-",checks)); dir.create(work)
        input <- file.path(work,"input.ancestrymap")
        output <- file.path(work,paste0("output.",format))
        if (!missing) writeBin(if(is.raw(text)) text else charToRaw(text),input)
        if (output_dir) dir.create(output) else writeBin(charToRaw("sentinel\n"),output)
        before <- sort(list.files(work))
        arguments <- list(input.file=input,output.file=output,force=TRUE)
        for (name in names(override)) arguments[name] <- override[name]
        job <- list(format=format,args=arguments,mode=mode,libs=.libPaths(),work=work)
        job_path <- file.path(root,paste0(checks,"-job.rds"))
        result_path <- file.path(root,paste0(checks,"-result.rds"))
        log_path <- file.path(root,paste0(checks,".log"))
        saveRDS(job,job_path)
        command <- list(command=file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe" else "Rscript"),
                        args=shQuote(c("--vanilla",child,job_path,result_path)),stdout=log_path,stderr=log_path)
        if ("timeout" %in% names(formals(system2))) command$timeout <- 15
        status <- do.call(system2,command)
        if (status!=0L || !file.exists(result_path))
            stop("Isolated ANCESTRYMAP job ",checks," (",format,") failed with status ",status,":\n",
                 paste(readLines(log_path,warn=FALSE),collapse="\n"))
        result <- readRDS(result_path)
        if (!is.null(pattern)) {
            stopifnot(identical(result$status,"error"),nzchar(result$message),grepl(pattern,result$message),
                      identical(sort(list.files(work)),before))
            if (!output_dir) stopifnot(identical(readBin(output,"raw",n=100),charToRaw("sentinel\n")))
        } else if (mode=="leak") {
            stopifnot(result$status=="success",result$count==400L,
                      is.na(result$before)||identical(result$before,result$after))
            cat("Repeated errors:",result$count,"; descriptors:",result$before,"->",result$after,"\n")
        } else {
            stopifnot(result$status=="success")
            value <- readChar(output,file.info(output)$size,useBytes=TRUE)
            if(.Platform$OS.type=="windows") value <- gsub("\r\n","\n",value,fixed=TRUE)
            stopifnot(identical(value,expected))
        }
        invisible(result)
    }
    for (format in c("geno","lfmm")) {
        for (bad in list(character(),NA_character_,c("first","second"),"",NULL)) {
            execute(format,override=list(output.file=bad),pattern="output.file")
            execute(format,override=list(input.file=bad),pattern="input.file")
        }
        execute(format,missing=TRUE,pattern="open ANCESTRYMAP input")
        execute(format,output_dir=TRUE,pattern="open ANCESTRYMAP output")
        for (bad in list(character(),""))
            execute(format,override=list(output.file=bad),mode="native",pattern="native filename")
        for (gt in c("nope","1junk","0.5","0x1","2147483648","-2147483649",strrep("9",100)))
            execute(format,paste0("a i ",gt,"\n"),pattern="integer genotype")
        for (gt in c("-2147483648","2147483647","+0001","0009")) {
            expected <- paste0(if(gt=="+0001")"1" else if(gt=="0009")"9" else gt,"\n")
            result <- execute(format,paste0("a i ",gt,"\n"),expected=expected)
            if (gt %in% c("-2147483648","2147483647"))
                stopifnot(sum(grepl("Warning: some genotypes",result$stdout))==1L)
        }
        execute(format,c(charToRaw("a i 1"),as.raw(0),charToRaw("hidden\n")),pattern="NUL")
        execute(format,"a i\n",pattern="missing genotype")
        execute(format,paste0("a ",strrep("I",20010)," 1\n"),expected="1\n")
        execute(format,paste0("a i 1 ",strrep("x",21000)," ignored annotation"),expected="1\n")
        execute(format,"a i 0\nb i 1",expected=if(format=="geno")"0\n1\n" else "0 1\n")
        for (bytes in c(510L,511L,512L)) {
            name <- paste0(strrep("\u00e9",255),if(bytes==511L)"a" else if(bytes==512L)"aa" else "")
            execute(format,charToRaw(enc2utf8(paste0(name," i 1\n"))),
                    pattern=if(bytes==512L)"maximum is 511 bytes" else NULL,expected="1\n")
        }
    }
    execute("geno",mode="leak")
    cat("ANCESTRYMAP isolated cleanup:",checks,"jobs passed\n")
})
