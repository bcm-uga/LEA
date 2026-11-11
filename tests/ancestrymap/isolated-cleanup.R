# Child-only entry point: errors/crashes must never kill the parent test process.
args <- commandArgs(TRUE)
job <- readRDS(args[1])
.libPaths(job$libs)
library(LEA)
setwd(job$work)
tryCatch({
    if (identical(job$mode, "leak")) {
        descriptor_dir <- if (dir.exists("/dev/fd")) "/dev/fd" else if (dir.exists("/proc/self/fd")) "/proc/self/fd" else NULL
        descriptors <- function() if (is.null(descriptor_dir)) NA_integer_ else length(list.files(descriptor_dir, all.files=TRUE, no..=TRUE))
        before <- descriptors()
        inputs <- c(paste0(strrep("a",512)," i 1\n"),
                    "a i 0\na j 1\nb i 2\nb j 9\nb k 0\nc i 1\n")
        count <- 0L
        for (input in inputs) {
            writeBin(charToRaw(input), "input.ancestrymap")
            for (i in seq_len(200)) {
                failure <- tryCatch(ancestrymap2geno("input.ancestrymap", "output.geno"), error=identity)
                stopifnot(inherits(failure,"error")); count <- count + 1L
            }
        }
        gc()
        after <- descriptors()
        stopifnot(is.na(before) || identical(before, after))
        result <- list(status="success", count=count, before=before, after=after)
    } else {
        log <- capture.output({
            if (identical(job$mode,"native")) {
                value <- .C(paste0("R_ancestrymap2",job$format),job$args$input.file,
                            job$args$output.file,N=0L,M=0L,PACKAGE="LEA")
            } else value <- do.call(if(job$format=="geno") ancestrymap2geno else ancestrymap2lfmm,job$args)
        })
        result <- list(status="success", value=value, stdout=log)
    }
}, error=function(e) result <<- list(status="error", message=conditionMessage(e)))
saveRDS(result,args[2])
