#create Class sNMF Project
setClass("snmfProject",
    slots = c(snmfProject.file = "character", projDir = "character", snmfDir = "character", 
        input.file = "character", runs = "list", K="integer",
        snmfClass.files = "vector", n="integer", L="integer",
        creationTime = "POSIXct")
)



# addRun

setGeneric("addRun.snmfProject", function(project="snmfProject", 
    run="snmfClass") attributes("snmfProject"));
setMethod("addRun.snmfProject", signature(project="snmfProject", 
    run="snmfClass"),
    function(project, run) {
        project@runs[[length(project@runs) + 1]] = run
        project@K = c(project@K, run@K)
        project@snmfClass.files = c(project@snmfClass.files, paste(run@directory, run@snmfClass.file, sep=""))

        return(project)
    }
)

# getRuns

setGeneric("getRuns.snmfProject", function(object, ...) 
    standardGeneric("getRuns.snmfProject"));
setMethod("getRuns.snmfProject", "snmfProject",
    function(object, K) {
        #check of K
        if (missing(K)) {
            K = object@K;
        } else if (!(all(K %in% object@K))) {
            stop("Unknown K!")
        } 
        K = unique(K)
        res = list()
        # check of the run number
        for (ku in K) {
            for (rep in which(object@K == ku))
                res[[length(res) + 1]] = object@runs[[rep]]
        }
        res
    }
)

# plot

# display lambda for a value of d, and a Manhattan plot for a value of K. 
setMethod("plot", "snmfProject",
    function(x, ...){
        s = summary(x);
        axe = NULL
        K = sort(unique(x@K))
        
        for (k in 1:length(K)) {
            tK = FALSE
            for (w in which(x@K == K[k])) {
                if (x@runs[[w]]@entropy)
                    tK = TRUE;
            }
            if (tK)
                axe = c(axe, K[k])        
        }

        plot(axe, s$crossEntropy[1,], ylab="Cross-entropy", 
                xlab="Number of ancestral populations",     ...)
    }
)

# G

setGeneric("G",  function(object, K, run) matrix)
setMethod("G", "snmfProject",
    function(object, K, run) {
        
        #check K
        if (missing(K)) {
            # if only one, that is the one
            if (length(unique(object@K)) == 1) {
                K = unique(object@K)
            } else {
                stop("Please choose a value of K among: ", 
                    paste(unique(object@K), collapse=" "))
            }
        } else {
            if (length(K) > 1) {
                stop("K is")  ## TODO 
            }
            K = test_integer("K", K, NULL)
            if (!(K %in% unique(object@K))) {
                stop(paste("No run exists for K = ", K,
                    ". Please, choose a value of K among: ", 
                    paste(unique(object@K), collapse=" "),sep=""))
            }
        }

        # check of run
        r = which(object@K == K)
        if (missing(run)) {
            if (length(r) > 1) {
                stop(paste(length(r)," runs have been performed for K =", K,
                    ".\n", "Please choose one with the paramater 'run'"))
            } else {
                run = 1;
            }
        } else { 
            run = test_integer('run', run, NULL)
            if (run > length(r)) {
                stop(paste("You chose run number ", run,". But only ", 
                length(r)," run(s) have been performed.", sep=""))
            }
        } 
            
        R = Gvalues(object@runs[[r[run]]], paste(object@projDir, object@snmfDir, sep = ""))

        return(R)
    }
)

# Q

setGeneric("Q",  function(object, K, run) matrix)
setMethod("Q", "snmfProject",
    function(object, K, run) {
        
        #check K
        if (missing(K)) {
            # if only one, that is the one
            if (length(unique(object@K)) == 1) {
                K = unique(object@K)
            } else {
                stop("Please, choose a value of K among: ", 
                    paste(unique(object@K), collapse=" "))
            }
        } else {
            K = test_integer('K', K, NULL)
            if (!(K %in% unique(object@K))) {
                stop(paste("No run exists for K = ", K,
                    ". Please, choose a value of K among: ", 
                    paste(unique(object@K), collapse=" "),sep=""))
            }
        }
        # check of run
        r = which(object@K == K)
        if (missing(run)) {
            if (length(r) > 1) {
                stop(paste(length(r)," runs have been performed for K =", K,
                    ".\n","Please choose one with the paramater 'run'", sep=""))
            } else {
                run = 1;
            }
        } else {
            run = test_integer('run', run, NULL)
            if (run > length(r)) {
                stop(paste("You chose run number ", run,". But only ", 
                length(r)," run(s) have been performed.", sep=""))
            }
        } 
            
        R = Qvalues(object@runs[[r[run]]], paste(object@projDir, object@snmfDir, sep = ""))

        return(R)
    }
)

# crossEntropy

setGeneric("cross.entropy", function(object, K, run) vector)
setMethod("cross.entropy", "snmfProject",
    function(object, K, run) {

        #check of K
        if (missing(K)) {
            # if only one, that is the one
            if (length(unique(object@K)) == 1) {
                K = unique(object@K)
            } else {
                stop(paste("Please, choose a value of K among: ", 
                    paste(unique(object@K), collapse=" "), sep=""))
            }
        } else if (!(K %in% unique(object@K))) {
            stop(paste("No run exists for K = ", K,
                ". Please, choose a value of K among: ", 
                paste(unique(object@K), collapse=" "),sep=""))
        }

        # check of run
        r = which(object@K == K)
        if (missing(run)) {
            run = 1:length(r);
        } 

        colnames = paste("K =", K)
        rownames = NULL
        res = NULL
        for (i in run) {
            if (i > length(r)) {
                stop(paste("You chose run number ", i,". But only ", 
                length(r)," run(s) have been performed.", sep=""))
            }
            if (object@runs[[r[i]]]@entropy) { 
                rownames = c(rownames, paste("run", i)) 
                res = c(res, object@runs[[r[i]]]@crossEntropy)
            }
        } 
            
        if (length(res) == 0) {
            cat(paste("The selected runs are without cross-entropy", 
                "criterion estimation!\n"))
            return(NULL);
        } else {
            return(matrix(res, dimnames = list(rownames, colnames)));
        }
    }
) 


# snmf-pvalues
setGeneric("snmf.pvalues", function(object, genomic.control, lambda, ploidy, entropy, fisher,  K, run) list)
setMethod("snmf.pvalues", "snmfProject",
          function(object, genomic.control, lambda, ploidy, entropy, fisher, K, run) {
           
             # check for fisher test
            if (missing(fisher)) {
              fisher <- TRUE
            }

           # check for genomic control
            if (missing(genomic.control)) {
              genomic.control <- TRUE
            }
            
            # check for lambda
            if (missing(lambda)) {
              lambda <-  1.0
            }
  

            # check for run
 
            if (missing(run)) {
              if (entropy) {
                ce <-  cross.entropy(object, K)
                run <-  which.min(ce) 
                } else {
                  stop("Run number is missing.") 
                }
            }
            
            if (length(run) > 1) {
              stop("Select a single run number.")
            }
            
            # check for ploidy
            if (missing(ploidy)) {
              stop("The ploidy parameter is missing without any default value.")
            }
            if (ploidy > 2) {
              stop("Choose option among ploidy = 1 or ploidy = 2.")
            }
            
            # estimate Fst            
            l <-  nrow(G(object, K, run))
            q <-  apply(Q(object, K, run), MARGIN = 2, mean)
            
            if (ploidy == 2) {
                G1.t <-  G(object, K, run)[seq(2,l,by = 3),]
                G2.t <-  G(object, K, run)[seq(3,l,by = 3),]
                freq <-  G1.t/2 + G2.t
              } else {
                freq <-  G(object, K, run)[seq(2,l,by = 2),]
              }
            
            H.s <-  apply(freq*(1-freq), MARGIN = 1, FUN = function(x) sum(q*x))
            P.t <-  apply(freq, MARGIN = 1, FUN = function(x) sum(q*x))
            H.t <-  P.t*(1-P.t)
            fst <-  1 - H.s/H.t
            
            
            # compute the z-scores
            
            n <-  nrow(Q(object, K, run))
            
            fst[fst<0] <-  0.000001
            
            zsq <- fst*(n-K)/(1-fst)
            
            if (genomic.control) {
              if (!fisher) {
              gif <- median(zsq)/qchisq(0.5, df = K-1) } else {
              gif <- median(zsq)/qf(0.5, df1 = K-1, df2 = n-K) }
              }
           else {
              gif <- lambda }
            
            if (!fisher) {
              snmf.pvalues <-  pchisq(zsq/gif, df = K-1, lower.tail = FALSE)}
            else {
              snmf.pvalues <-  pf(zsq/gif, df1 = K-1, df2 = n-K, lower.tail = FALSE)
            }
            
            return(list(pvalues = snmf.pvalues, GIF = gif))
          }
)

#impute

setGeneric("impute",  function(object, input.file, method, K, run) NULL)
setMethod("impute", "snmfProject",
          function(object, input.file, method, K, run) {
            
            
            # check input file
            
            
            if (missing(input.file)){
              stop("Missing input file: geno or lfmm format\n")} else {
                input.file = test_character("input.file", input.file, NULL)
                # check extension and convert if necessary 
                input.file = test_input_file(input.file, "lfmm") 
              }
            if (!file.exists(input.file)) stop("Input file not found.\n")
            
            
            #check  K
            if (missing(K)) {
              # if only one, that is the one
              if (length(unique(object@K)) == 1) {
                K = unique(object@K)
              } else {
                stop("Please, choose a value of K among: ", 
                     paste(unique(object@K), collapse = " "))
              }
            } else {
              if (length(K) > 1) {
                stop("K is")  ## TODO 
              }
              K = test_integer("K", K, NULL)
              cat(paste("Missing genotype imputation for K =",K,"\n"))
              
              if (!(K %in% unique(object@K))) {
                stop(paste("No run exists for K = ", K,
                           ". Please, choose a value of K among: ", 
                           paste(unique(object@K), collapse=" "),sep=""))
              }
            }
            
            # check of run
            r = which(object@K == K)
            if (missing(run)) {
              if (length(r) > 1) {
                stop(paste(length(r)," runs have been performed for K =", K,
                           ".\n", "Please choose one run with the parameter 'run'"))
              } else {
                run = 1;
              }
            } else { 
              run = test_integer('run', run, NULL)
              cat(paste("Missing genotype imputation for run =",run,"\n"))
              
              
              if (run > length(r)) {
                stop(paste("You chose run number ", run,"but only ", 
                           length(r),"run(s) were performed.", sep=""))
              }
            } 
            
            
            # check of method
            if (missing(method)) {
              method = "mode" 
              } else { 
                method = test_character('method', method, NULL)
                if ( sum(method == c("random","mode")) < 1 ) 
                  {stop("method must be 'mode' or 'random'.")}
              } 
            
            
            dat <- read.lfmm(input.file)
            if (!(9 %in% dat)) {stop("Missing data must be encoded as 9.")} 
            QQ <- Q(object, K, run)
            GG <- t(G(object, K, run))
            nploidy <- ncol(GG)/object@L 
            
            X <-  QQ %*% GG
            
            for (i in 1:nrow(dat)){
              mat <- matrix(X[i,], nrow = nploidy, byrow = FALSE)
              
              if (method == 'mode') {
                wm <- apply(mat, 2, which.max) 
                } else { 
                  wm <- apply(mat, 2, FUN = function(x) sample(1:nploidy, 1, prob = x))
                }
              
              dat[i, dat[i,]==9] <- wm[dat[i,] == 9] - 1  
              
            }
            out.file  <- paste(input.file,"_imputed.lfmm", sep ="")
            
            write.lfmm(R = dat, output.file = out.file )
            
            rm(X)
            rm(dat)
            cat("Results are written in the file: ",
                out.file,"\n" )
          }
)



#barchart
setGeneric("barchart",  function(object, K, run, sort.by.Q = TRUE, lab = FALSE,...) list)
setMethod("barchart", "snmfProject",
          function(object, K, run, sort.by.Q = TRUE, lab = FALSE,...){
            
            #check  K
            if (missing(K)) {
              # if only one, that is the one
              if (length(unique(object@K)) == 1) {
                K = unique(object@K)
              } else {
                stop("Please choose a value of K among: ", 
                     paste(unique(object@K), collapse=" "))
              }
            } else {
              if (length(K) > 1) {
                stop("K is")  ## TODO 
              }
              K = test_integer("K", K, NULL)
              
              if (!(K %in% unique(object@K))) {
                stop(paste("No run exists for K = ", K,
                           ". Please, choose a value of K among: ", 
                           paste(unique(object@K), collapse=" "),sep=""))
              }
            }
            
            # check of run
            r = which(object@K == K)
            if (missing(run)) {
              if (length(r) > 1) {
                stop(paste(length(r)," runs have been performed for K =", K,
                           ".\n", "Please choose one run with the parameter 'run'"))
              } else { run = 1 }
            } else { 
              run = test_integer('run', run, NULL)

              if (run > length(r)) {
                stop(paste("You chose run number ", run,". But only ", 
                           length(r)," run(s) have been performed.", sep=""))
              }
            } 
            
            QQ <- Q(object, K, run)
            
            if (sort.by.Q) {
              gr <- apply(QQ, MARGIN = 1, which.max)
              gm <-  max(gr)
              gr.o <-  order(sapply(1:gm, FUN = function(g) mean(QQ[,g])))
              gr <- sapply(gr, FUN = function(i) gr.o[i])
              or <-  order(gr)
              Qm <-  t(QQ[or,])
              class(Qm) <-  "matrix"
              graphics::barplot(height = Qm, ...)
              return(list(order = or))
            } else {
              Qm <-  t(QQ)
              class(Qm) <-  "matrix"
              graphics::barplot(height = Qm, ...)
              return(list(order = 1:nrow(QQ)))
            }
          }
)







# show

setMethod("show", "snmfProject",
    function(object) {
        cat("snmf Project\n\n")
        cat("snmfProject file:                ", object@snmfProject.file, "\n")
        cat("project directory:               ", object@projDir, "\n")
        cat("snmf results directory:          ", object@snmfDir, "\n")
        cat("date of creation:                ", object@creationTime, "\n")
        cat("input file:                      ", object@input.file, "\n")
        cat("number of individuals:           ", object@n, "\n")
        cat("number of loci:                  ", object@L, "\n")
        cat("number of ancestral populations: ", object@K, "\n")
        if (length(object@runs)) {
            for (i in 1:length(object@runs)) {
                cat("\n")
                cat("***** run *****\n");
                show(object@runs[[i]])
            }
        }
    }
)

# summary

setGeneric("summary", function(object) NULL)
setMethod("summary", "snmfProject",
    function(object) {
        K = sort(unique(object@K))
        rownames=c("with cross-entropy", "without cross-entropy", "total")
        colnames=paste("K =", K)
        rep = matrix(NA, ncol=length(K), nrow=3, 
            dimnames= list(rownames, colnames))
        for (k in 1:length(K)) {
            rep[3,k] = length(which(object@K == K[k]))
            rep[1,k] = 0;
            for (w in which(object@K == K[k])) {
                if (object@runs[[w]]@entropy)
                    rep[1,k] = rep[1,k] + 1
            }
            rep[2,k] = rep[3,k] - rep[1,k]
        }
        
        rownames = c("min", "mean", "max");
        ce = matrix(NA, ncol=length(K), nrow=3, 
            dimnames= list(rownames, colnames)); 
        for (k in 1:length(K)) {
            ceK = cross.entropy(object, K[k])
            if (!is.null(ceK)) { 
                ce[1,k] = min(ceK);
                ce[2,k] = mean(ceK);
                ce[3,k] = max(ceK);
            } else {
                ce[1,k] = NA
                ce[1,k] = NA
                ce[1,k] = NA
            }
        }
    list(repetitions=rep, crossEntropy=ce)
    }
)

# load

setGeneric("load.snmfProject", function(file="character") 
    attributes("snmfProject"))
setMethod("load.snmfProject", "character",
    function(file) {
        res = dget(file);
        if (length(res@snmfClass.files) > 0) {
            for (r in 1:length(res@snmfClass.files)) {
                res@runs[[r]] = load.snmfClass(paste(res@projDir, res@snmfDir, 
                    res@snmfClass.files[r], sep = ""))
            }
        }
        return(res);
    }
)

# save

setGeneric("save.snmfProject", function(object) character)
setMethod("save.snmfProject", signature(object="snmfProject"),
    function(object) {
        file = object@snmfProject.file;
        if (length(object@runs) > 0) {
            for (r in 1:length(object@runs)) {
                save.snmfClass(object@runs[[r]], paste(object@projDir, object@snmfDir,
                    object@snmfClass.files[r], sep = ""))
            }
        }
        object@runs = list()
        pathFile = paste(object@projDir, file, sep = "")
        dput(object, pathFile) 
        cat("The project is saved into :\n",currentDir(pathFile),"\n\n");
        cat("To load the project, use:\n project = load.snmfProject(\"",
            currentDir(pathFile),"\")\n\n",sep="");
        cat("To remove the project, use:\n remove.snmfProject(\"",
            currentDir(pathFile),"\")\n\n",sep="");
   
        pathFile
    }
)

# remove

setGeneric("remove.snmfProject", function(file="character") NULL)
setMethod("remove.snmfProject", "character",
    function(file) {
        # file
        test_character("file", file, NULL);

        res = dget(file);
        unlink(paste(res@projDir, res@snmfDir, sep = ""), recursive = TRUE)
        file.remove(file)
    }
)

# export

setGeneric("export.snmfProject", function(file, force) character)
setMethod("export.snmfProject", "character",
    function(file, force) {
        # file 
        test_character("file", file, NULL);
        # force
        force = test_logical("force", force, FALSE)

        object = load.snmfProject(file)

        pathFile = paste(object@projDir, object@snmfProject.file, sep = "")
        zipFile = paste(setExtension(pathFile, ""), "_snmfProject.zip", sep = "")

        if (force == FALSE && file.exists(zipFile)) {
            stop("An export with the same name already exists.\n",
                 "If you want to overwrite it, please use the ",
                 "option 'force == TRUE'\n\n") 
        } else {
            if (file.exists(zipFile))
                file.remove(zipFile)
            curDir = getwd()
            setwd(object@projDir)
            zip(zipFile, c(object@snmfProject.file,
                paste(object@snmfDir, sep = ""), object@input.file))
            setwd(curDir)
            cat("An export of the snmf project hase been created: ", currentDir(zipFile), "\n\n") 
        } 
        
    }
)

# import

setGeneric("import.snmfProject", function(zipFile, directory, force) attributes("snmfProject"))
setMethod("import.snmfProject", "character",
    function(zipFile, directory, force) {
        #file 
        zipFile = test_character("zipfile", zipFile, NULL);
        # directory
        directory = test_character("directory", directory, getwd())
        # force
        force = test_logical("force", force, FALSE)

        # check that no file exists
        tmp = basename(zipFile)
        file = paste(normalizePath(directory), "/", substr(tmp, 1, 
            nchar(tmp) - 16), ".snmfProject", sep = "")
        if (!force && file.exists(file)) {
            stop("A snmf project with same name exists in the directory: \n",
                 directory, ".\n",
                 "If you want to overwrite it, please use the ",
                 "option 'force == TRUE'\n\n") 

        }
        # unzip
        unzip(zipFile, exdir = directory)
        cat("The project has been imported into directory, ", directory, "\n\n") 
        # modify the project infos
        res = dget(file);
        res@projDir = paste(normalizePath(directory), "/", sep = "")
        res@creationTime = Sys.time()
        dput(res, file)
        # load the project
        load.snmfProject(file)
    } 
)

# combine

setGeneric("combine.snmfProject", function(combined.file, file.to.combine, force) attributes("snmfProject"))
setMethod("combine.snmfProject", signature(combined.file = "character", file.to.combine = "character"),
    function(combined.file, file.to.combine, force) {
        # file 
        test_character("combined.file", combined.file, NULL);
        # file 
        test_character("file.to.combine", file.to.combine, NULL);
        # force
        force = test_logical("force", force, FALSE)

        to.combine = load.snmfProject(file.to.combine)
        combined = load.snmfProject(combined.file)

        # check that the projects are compatible
        # check n
        if (to.combine@n != combined@n) {
            stop("The number of individuals (",combined@n," and ",to.combine@n,") are\n",
            "different in the two projects. The two following projects cannot be\n",
            "combined:", combined.file, "\n", file.to.combine, "\n\n")
        # check L    
        } else if (to.combine@L != combined@L) {
            stop("The number of loci (",combined@L," and ",to.combine@L,") are\n", 
            "different in the two projects. The two following projects cannot be\n",
            "combined:", combined.file, "\n", file.to.combine, "\n\n")
        } 
        # check input file name
        if (!force && basename(to.combine@input.file) != basename(combined@input.file)) { 
            stop("The names of the input file of the two projects are different \n",
                 ": ", basename(to.combine@input.file), ", ", basename(combined@input.file), "\n",
                 "Be cautious that only projects with the same input file can be\n",
                 "combined. If you still want to combine them, please use the ",
                 "option 'force == TRUE'\n\n") 
        }

        # let's go  
        if (length(to.combine@runs) > 0) { 
            for (r in 1:length(to.combine@runs)) {
                to.combine.run = to.combine@runs[[r]]
                combined.run = to.combine@runs[[r]]
                k = combined.run@K
                re = as.integer(length(which(combined@K == k)) + 1)
                # combine  
                # create file directory
                tmp  = basename(setExtension(basename(combined@input.file), ""))
                combined.run@directory = paste("K", k, "/run", re, "/", sep="") 
                dir.create(paste(combined@projDir, combined@snmfDir, combined.run@directory, 
                    sep = ""), showWarnings = FALSE, recursive = TRUE)
                # snmfClass.file
                combined.run@snmfClass.file = paste(tmp, "_r", re ,".",k, 
                    ".snmfClass", sep = "")
                # run
                combined.run@run = re
                # Q
                combined.run@Q.output.file = paste(tmp, "_r", re ,".",k, ".Q", sep="")
                file.copy(paste(to.combine@projDir, to.combine@snmfDir, to.combine.run@directory, 
                    to.combine.run@Q.output.file, sep = ""), 
                    paste(combined@projDir, combined@snmfDir, combined.run@directory, 
                    combined.run@Q.output.file, sep = "")) 
                # G
                combined.run@G.output.file = paste(tmp, "_r", re ,".",k, ".G", sep="")
                file.copy(paste(to.combine@projDir, to.combine@snmfDir, to.combine.run@directory, 
                    to.combine.run@G.output.file, sep = ""), 
                    paste(combined@projDir, combined@snmfDir, combined.run@directory, 
                    combined.run@G.output.file, sep = "")) 
                # save run
                save.snmfClass(combined.run, paste(combined@projDir, combined@snmfDir,
                    combined.run@directory, combined.run@snmfClass.file, sep = ""))
                # save project
                combined = addRun.snmfProject(combined, combined.run)
            }
        }
        save.snmfProject(combined)
    }
)
