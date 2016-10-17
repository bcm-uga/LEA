# create Class LFMM Project
setClass("lfmmProject",
    slots = c(lfmmProject.file = "character", projDir="character", lfmmDir = "character",
        input.file = "character", environment.file = "character", 
        runs = "list", K="integer", d="integer", all = "logical",
        lfmmClass.files = "vector", n="integer", L = "integer",
        D = "integer", creationTime = "POSIXct")
)

# addRun

setGeneric("addRun.lfmmProject", function(project="lfmmProject", 
    run="lfmmClass") attributes("lfmmProject"));
setMethod("addRun.lfmmProject", signature(project="lfmmProject", 
    run="lfmmClass"), function(project, run) {
        project@runs[[length(project@runs) + 1]] = run
        project@K = c(project@K, run@K)
        project@d = c(project@d, run@d)
        project@all = c(project@all, run@all)
        project@lfmmClass.files = c(project@lfmmClass.files, 
            paste(run@directory, run@lfmmClass.file, sep = ""));

        return(project)
    }
)

# adjusted-pvalues

setGeneric("adjusted.pvalues", function(object, genomic.control, lambda, K, d, all, run) vector)
setMethod("adjusted.pvalues", "lfmmProject",
    function(object, genomic.control, lambda, K, d, all, run) {

        # check for genomic control
        if (missing(genomic.control)) {
            genomic.control = TRUE
        }

        # check for lambda
        if (missing(lambda)) {
            lambda = 1.0
        }

        # calculate the z-scores
        res = z.scores(object, K, d, all, run)

        zs.median <- apply(res, MARGIN = 1, median)
    
        if (genomic.control)
            gif <- median(zs.median^2)/qchisq(0.5, df=1)
        else
            gif <- lambda

        adjusted.pvalues = pchisq(zs.median^2/gif, df = 1, lower.tail = FALSE)

        return(list(p.values = adjusted.pvalues, genomic.inflation.factor = gif))
    }
)


# listMethods

#listMethods_lfmmProject <- function()
#{
#}

# listSlots

#listSlots_lfmmProject <- function()
#{
#c(    "lfmmProject_file", "input_file", "runs", "K", "d", "all", 
# "lfmmClass_files")
#}

#.DollarNames.lfmmProject <- function(x, pattern) c(listSlots_lfmmProject(),
#    listMethods_lfmmProject())

# $

#setMethod("$", "lfmmProject",
#           function(x, name) {
#             if (!(name %in% listMethods_lfmmProject() || name %in% 
#               listSlots_lfmmProject())) {
#               stop("no $ method for object without attributes")
#         } else if (name %in% listMethods_lfmmProject()) {
#        do.call(name, list(x));
#         } else {
#        slot(x, name)
#         }
#           }
#)

# getRuns

setGeneric("getRuns.lfmmProject", function(object, k, d, all) 
    standardGeneric("getRuns.lfmmProject"));
setMethod("getRuns.lfmmProject", "lfmmProject",
    function(object, k, d, all) {
        
        # check of k
        if (missing(k)) {
            k = object@K
        } else if (!(all(k %in% object@K))) {
            stop("Unknown k!")
        }
        k = unique(k)
        runk = NULL
        for (ku in k) 
            runk = c(runk, which(object@K == ku))
        
        # check d
        if (missing(d)) {
            d = object@d
        } else if (!(all(d %in% object@d))) {
            stop("Unknown d!")
        }
        d = unique(d)
        rund = NULL
        for (du in d) 
            rund = c(rund, which(object@d == du))

        # check all    
        if (missing(all)) {
            all = c(TRUE, FALSE)
        }
        runall = NULL
        for (allu in all) 
            runall = c(runall, which(object@all == allu))
        
        rep = intersect(runk, rund)
        rep = intersect(rep, runall)

        res = list()
        for (r in rep) {
        # check of the run number
        res[[length(res) + 1]] = object@runs[[r]]
        }
    }
)

# show

setMethod("show", "lfmmProject",
    function(object) {
        cat("lfmm Project\n\n")
        cat("lfmmProject file:                     ", object@lfmmProject.file,
            "\n")
        cat("project directory:                    ", object@projDir, "\n")
        cat("lfmm result directory:                ", object@lfmmDir, "\n")
        cat("date of creation:                     ", object@creationTime,
            "\n")
        cat("input file:                           ", object@input.file, "\n")
        cat("variable file:                        ", object@environment.file,
            "\n")
        cat("number of individuals:                ", object@n, "\n")
        cat("number of loci:                       ", object@L, "\n")
        cat("number of environmental variables:    ", object@D, "\n")
        cat("run for the latent factors:           ", object@K, "\n")
        cat("run for the environmental variables:  ", object@d, "\n")
        cat("variables separately or together:     ", object@all, "\n")
        if (length(object@runs)) {
            for (i in 1:length(object@runs)) {
                cat("\n")
                cat("****** run *******\n");
                show(object@runs[[i]])
            }
        }
    }
)

# summary

setGeneric("summary", function(object) NULL)
setMethod("summary", "lfmmProject",
    function(object) {
        separately = which(!object@all)
        together = which(object@all)
        d_separately = unique(object@d[separately])
        d_together = unique(object@d[together])

        # dimnames
        if (length(d_separately) > 0) 
            rownames_sep = paste("separately d =", d_separately)
        else 
            rownames_sep = NULL
        if (length(d_together) > 0) 
            rownames_tog = paste("together d =", d_together)
        else 
            rownames_tog = NULL

        K = sort(unique(object@K))
        rownames=c(rownames_sep, rownames_tog, "total")
        colnames=paste("K =", K)
        rep = matrix(0, ncol=length(colnames), nrow=length(rownames), 
            dimnames= list(rownames, colnames))
        # fill repetitions
        for (k in 1:length(K)) {
            total = 0;
            r = which(object@K == K[k])
            if (length(d_separately) > 0) { 
                for(d in 1:length(d_separately)) {
                    rep[d,k] = length(intersect(r, 
                        which(object@d[separately] == d_separately[d])))
                    total = total + rep[d,k]
                }
            }
            if (length(d_together) > 0) { 
                for(d in 1:length(d_together)) {
                    rep[length(d_separately) + d,k] = length(intersect(r, 
                        which(object@d[together] == d_together[d])))
                    total = total + rep[d,k]
                }
            }
            rep[length(d_separately)+length(d_together) + 1,k] = total
        }
        
#        # fill the inflation factor
#        rownames = c(rownames_sep, rownames_tog)
#        la = matrix(NA, ncol=length(colnames), nrow=length(rownames), 
#            dimnames= list(rownames, colnames))
#        # fill repetitions
#        for (k in 1:length(K)) {
#            r = which(object@K == K[k])
#            if (length(d_separately) > 0) { 
#                for(d in 1:length(d_separately)) {
#                    inter = intersect(r,
#                     which(object@d[separately] == d_separately[d]));
#                    l = NULL
#                    for (it in inter)
#                        l = c(l, object@runs[[it]]@inflationFactor)
#                    if (!is.null(l))
#                        la[d,k] = mean(l)
#                }
#            }
#            if (length(d_together) > 0) { 
#                for(d in 1:length(d_together)) {
#                    inter = intersect(r,
#                       which(object@d[together] == d_together[d]));
#                    l = NULL
#                    for (it in inter)
#                        l = c(l, object@runs[[it]]@inflationFactor)
#                    if (!is.null(l))
#                        la[length(d_separately) + d,k] = mean(l)
#                }
#            }
#        }
        list(repetitions= rep) #, inflationFactor= la)
    }
)


# cFR
#setGeneric("cFDR", function(object, percentage, k, d, all, cutoff) vector)
#setMethod("cFDR", c("lfmmProject", "numeric"),
#    function(object, percentage, k, d, all, cutoff) {
#
#        # check of k
#        if (missing(k) && length(unique(object@K)) == 1) {
#            k = object@K
#        } else if (missing(k) || !k %in% object@K) {
#            stop("Unknown k!")
#        }
#        runk = which(object@K == k)
#        
#        # check d
#        if (missing(d) && length(unique(object@d)) == 1) {
#            d = object@d
#        } else if (missing(d) || !(all(d %in% object@d))) {
#            stop("Unknown d!")
#        }
#        rund = which(object@d == d)
#
#        # check all    
#        if (missing(all) && length(unique(object@all)) == 1) {
#            all = object@all
#        } else if (missing(all) || !(all %in% object@all)) {
#            stop("Unknown all!")
#        }
#        runall = which(object@all == all)
#        
#        # check cutoff
#        if (missing(cutoff)) {
#            cutoff = 0;
#        }
#
#        rep = intersect(runk, rund)
#        rep = intersect(rep, runall)
#
#        lst = NULL
#        for (r in rep) {
#            lst = c(lst, cFDR(object@runs[[r]],percentage))
#        }
#        tab1 = cbind(sort(unique(lst)),table(lst))
#        res = tab1[tab1[,2]>cutoff,1]
#
#        return(res);
#    } 
#)

# plot
#setMethod("plot", "lfmmProject",
#          function(x, ...){
#        s = summary(x)  
#        par(mfrow=c(1,dim(s@inflationFactor)[1]))
#        for (i in 1:dim(s$inflationFactor[1])) {
#           plot(s$inflationFactor[i,], ylab="Cross-Entropy", 
#           xlab="Number of ancestral populations", ...)
#        }
#          }
#)

# mlog10pvalues

setGeneric("mlog10p.values", function(object, K, d, all, run) vector)
setMethod("mlog10p.values", "lfmmProject",
    function(object, K, d, all, run) {
        # check of k
        if (missing(K)) {
            if (length(unique(object@K)) == 1) {
                K = unique(object@K)
            } else {
                stop(paste("Please, choose a value of K among: ", 
                    paste(unique(object@K), collapse=" "), sep=""))
            }
        } else if (!(K %in% object@K)) {
            stop(paste("No run exists for K = ",K,
                ". Please, choose a value of K among:", paste(unique(object@K),
                collapse=" "), sep = ""))
        }
        runk = which(object@K == K)

        # check d
        if (missing(d)) {
            if(length(unique(object@d[runk])) == 1) {
                d = unique(object@d[runk])
            } else {
                stop(paste("Please, choose a value of d among: ", 
                paste(unique(object@d[runk]), collapse=" "), sep=""))
            }
        } else if (!(d %in% object@d[runk])) {
            stop(paste("Please, choose a value of d among: ", 
                paste(unique(object@d[runk]), collapse=" "), sep=""))
        }
        rund = intersect(which(object@d == d), runk)

        # check all  
        if(missing(all)) {
            if(length(unique(object@all[rund])) == 1) {
                all = unique(object@all[rund])
            } else {
                stop(paste("Please, choose a value of the all parameter among:",
                paste(unique(object@all[rund]), collapse=" "), sep=""))
            }
        }

        # check run 
        r = intersect(which(object@all == all), rund)
        if (missing(run))
            run = 1:length(r)

        colnames = paste("run", run)
        rownames = NULL
        res = matrix(NA, ncol=length(run), nrow = object@L, 
            dimnames = list(rownames, colnames))
        for (i in 1:length(run)) {
            if (run[i] > length(r)) {
                stop(paste("You chose run number ", run[i],". But only ",
                    length(r)," run(s) have been performed.", sep=""))
            }
            res[,i] = mlog10pvalues(object@runs[[r[run[i]]]], 
                paste(object@projDir, object@lfmmDir, sep = ""))        
        }

        return(res)
    }
)

# pvalues

setGeneric("p.values", function(object, K, d, all, run) vector)
setMethod("p.values", "lfmmProject",
    function(object, K, d, all, run) {
        # check of k
        if (missing(K)) {
            if (length(unique(object@K)) == 1) {
                K = unique(object@K)
            } else {
                stop(paste("Please, choose a value of K among: ", 
                    paste(unique(object@K), collapse=" "), sep=""))
            }
        } else if (!(K %in% object@K)) {
            stop(paste("No run exists for K = ",K,
                ". Please, choose a value of K among:", 
                paste(unique(object@K), collapse=" "), sep = ""))
        }
        runk = which(object@K == K)

        # check d
        if (missing(d)) {
            if(length(unique(object@d[runk])) == 1) {
                d = unique(object@d[runk])
            } else {
                stop(paste("Please, choose a value of d among: ", 
                    paste(unique(object@d[runk]), collapse=" "), sep=""))
            }
        } else if (!(d %in% object@d[runk])) {
            stop(paste("Please, choose a value of d among: ", 
                paste(unique(object@d[runk]), 
                collapse=" "), sep=""))
                }
        rund = intersect(which(object@d == d), runk)

        # check all  
        if(missing(all)) {
            if(length(unique(object@all[rund])) == 1) {
                all = unique(object@all[rund])
            } else {
                stop(paste("Please, choose a value of the all parameter among:",
                    paste(unique(object@all[rund]), collapse=" "), sep=""))
            }
        }

        # check run 
        r = intersect(which(object@all == all), rund)
        if (missing(run))
            run = 1:length(r)

        colnames = paste("run", run)
        rownames = NULL
        res = matrix(NA, ncol=length(run), nrow = object@L, 
            dimnames = list(rownames, colnames))
        for (i in 1:length(run)) {
            if (run[i] > length(r)) {
                stop(paste("You chose run number ", run[i],". But only ",
                    length(r)," run(s) have been performed.", sep=""))
            }
            res[,i] = pvalues(object@runs[[r[run[i]]]], 
                paste(object@projDir, object@lfmmDir, sep = ""))        
        }

        return(res)
    }
)

# zscores

setGeneric("z.scores", function(object, K, d, all, run) vector)
setMethod("z.scores", "lfmmProject",
    function(object, K, d, all, run) {
        # check of k
        if (missing(K)) {
            if (length(unique(object@K)) == 1) {
                K = unique(object@K)
            } else {
                stop(paste("Please, choose a value of K among: ", 
                    paste(unique(object@K), collapse=" "), sep=""))
            }
        } else if (!(K %in% object@K)) {
            stop(paste("No run exists for K = ",K,
                ". Please, choose a value of K among:", 
                paste(unique(object@K), collapse=" "), sep = ""))
        }
        runk = which(object@K == K)

        # check d
        if (missing(d)) {
            if(length(unique(object@d[runk])) == 1) {
                d = unique(object@d[runk])
            } else {
                stop(paste("Please, choose a value of d among: ", 
                    paste(unique(object@d[runk]), collapse=" "), sep=""))
            }
        } else if (!(d %in% object@d[runk])) {
            stop(paste("Please, choose a value of d among: ", 
                paste(unique(object@d[runk]), collapse=" "), sep=""))
        }
        rund = intersect(which(object@d == d), runk)
    
        # check all  
        if(missing(all)) {
            if(length(unique(object@all[rund])) == 1) {
                all = unique(object@all[rund])
            } else {
                stop(paste("Please, choose a value of the all parameter among:",
                paste(unique(object@all[rund]), collapse=" "), sep=""))
            }
        }

        # check run 
        r = intersect(which(object@all == all), rund)
        if (missing(run))
            run = 1:length(r)

        colnames = paste("run", run)
        rownames = NULL
        res = matrix(NA, ncol=length(run), nrow = object@L, 
            dimnames = list(rownames, colnames))
        for (i in 1:length(run)) {
            if (run[i] > length(r)) {
                stop(paste("You chose run number ", run[i],". But only ",
                    length(r)," run(s) have been performed.", sep=""))
            }
            res[,i] = zscores(object@runs[[r[run[i]]]], 
                paste(object@projDir, object@lfmmDir, sep = ""))        
        }

        return(res)
    }
)

# load

setGeneric("load.lfmmProject", function(file="character") 
    attributes("lfmmProject"))
setMethod("load.lfmmProject", "character",
    function(file) {
        # file 
        test_character("file", file, NULL);

        # load the project
        res = dget(file);
        if (length(res@lfmmClass.files) > 0) {
            for (r in 1:length(res@lfmmClass.files)) {
                res@runs[[r]] = load.lfmmClass(paste(res@projDir, res@lfmmDir,
                    res@lfmmClass.files[r], sep = ""))
            }
        }
        return(res);
    }
)

# save
setGeneric("save.lfmmProject", function(object="lfmmProject") character)
setMethod("save.lfmmProject", signature(object="lfmmProject"),
    function(object) {
        file = object@lfmmProject.file;
        if (length(object@runs) > 0) {
            for (r in 1:length(object@runs)) {
                save.lfmmClass(object@runs[[r]], paste(object@projDir, 
                    object@lfmmDir, object@lfmmClass.files[r], sep = ""))
            }
        }
        object@runs = list()
        pathFile = paste(object@projDir, file, sep = "")
        dput(object, pathFile) 
        cat("The project is saved into :\n",currentDir(pathFile),"\n\n");
        cat("To load the project, use:\n project = load.lfmmProject(\"",
            currentDir(pathFile),"\")\n\n",sep="");
        cat("To remove the project, use:\n remove.lfmmProject(\"",
            currentDir(pathFile), "\")\n\n", sep="")

        pathFile
    }
)

# remove

setGeneric("remove.lfmmProject", function(file="character") NULL)
setMethod("remove.lfmmProject", "character",
    function(file) {
        # file
        test_character("file", file, NULL);

        # load and remove
        res = dget(file);
        unlink(paste(res@projDir, res@lfmmDir, sep = ""), recursive = TRUE)
        file.remove(file)
    }
)

# export

setGeneric("export.lfmmProject", function(file, force) character)
setMethod("export.lfmmProject", "character",
    function(file, force) {
        # file 
        test_character("file", file, NULL);
        # entropy
        force = test_logical("force", force, FALSE)

        object = load.lfmmProject(file)

        pathFile = paste(object@projDir, object@lfmmProject.file, sep = "")
        zipFile = paste(setExtension(pathFile, ""), "_lfmmProject.zip", sep = "")

        if (force == FALSE && file.exists(zipFile)) {
            stop("An export with the same name already exists.\n",
                 "If you want to overwrite it, please use the ",
                 "option 'force = TRUE'\n\n") 
        } else {
            if (file.exists(zipFile))
                file.remove(zipFile)
            curDir = getwd()
            setwd(object@projDir)
            zip(zipFile, c(object@lfmmProject.file,
                paste(object@lfmmDir, sep = ""), object@input.file, object@environment.file))
            setwd(curDir)
            cat("An export of the lfmm project hase been created: ", currentDir(zipFile), "\n\n") 
        } 
        
    }
)

# import

setGeneric("import.lfmmProject", function(zipFile, directory, force) attributes("lfmmProject"))
setMethod("import.lfmmProject", "character",
    function(zipFile, directory, force) {
        # file 
        zipFile = test_character("zipfile", zipFile, NULL);
        # directory
        directory = test_character("directory", directory, getwd())
        # force
        force = test_logical("force", force, FALSE)

        # check that no file exists
        tmp = basename(zipFile)
        file = paste(normalizePath(directory), "/", substr(tmp, 1, 
            nchar(tmp) - 16), ".lfmmProject", sep = "")
        if (!force && file.exists(file)) {
            stop("A lfmm project with same name exists in the directory: \n",
                 directory, ".\n",
                 "If you want to overwrite it, please use the ",
                 "option 'force == TRUE'\n\n") 
        }
        # unzip
        unzip(zipFile, exdir = directory)
        cat("The project has been imported into directory, ", directory, "\n\n") 
        # modify the project directory
        res = dget(file);
        res@projDir = paste(normalizePath(directory), "/", sep = "")
        res@creationTime = Sys.time()
        dput(res, file)
        # load the project
        load.lfmmProject(file)
    } 
)

# combine

setGeneric("combine.lfmmProject", function(combined.file, file.to.combine, force) attributes("lfmmProject"))
setMethod("combine.lfmmProject", signature(combined.file = "character", file.to.combine = "character"),
    function(combined.file, file.to.combine, force) {
        # file 
        test_character("combined.file", combined.file, NULL);
        # file 
        test_character("file.to.combine", file.to.combine, NULL);
        # force
        force = test_logical("force", force, FALSE)

        to.combine = load.lfmmProject(file.to.combine)
        combined = load.lfmmProject(combined.file)
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
        # check D 
        } else if (to.combine@D != combined@D) {
            stop("The number of environmental variables (",combined@D," and ",to.combine@D,") are\n", 
            "different in the two projects. The two following projects cannot be\n",
            "combined:", combined.file, "\n", file.to.combine, "\n\n")
        } 
        # check input file name
        if (!force && basename(to.combine@input.file) != basename(combined@input.file)) { 
            stop("The names of the input file of the two projects are different:\n", 
                 basename(to.combine@input.file), "\n, ", basename(combined@input.file), "\n",
                 "Be cautious that only projects with the same input file can be\n",
                 "combined. If you still want to combine them, please use the ",
                 "option 'force == TRUE'\n\n") 
        }
        # check environment file
        to.combine.C = read.env(paste(to.combine@projDir, to.combine@environment.file, sep = ""))
        combined.C = read.env(paste(combined@projDir, to.combine@environment.file, sep = ""))
        diff = sum(abs(to.combine.C - combined.C)/combined.C)/length(combined.C)
        if (!force && diff > 0.01) {
            stop("The two environment files seem to be different:\n",
                 basename(to.combine@input.file), "\n, ", basename(combined@input.file), "\n",
                 "Be cautious that only projects with the same input file can be\n",
                 "combined. If you still want to combine them, please use the ",
                 "option 'force == TRUE'\n\n") 
        }
        # let's go  
        if (length(to.combine@runs) > 0) { 

            list.re = NULL
            for (r in 1:length(to.combine@runs)) {
                list.re = c(list.re, length(combined@K) + to.combine@runs[[r]]@run)
            }

            for (r in 1:length(to.combine@runs)) {
                to.combine.run = to.combine@runs[[r]]
                combined.run = to.combine@runs[[r]]
                k = combined.run@K
                re = list.re[r]
                # combine  
                # create file directory
                tmp  = basename(setExtension(basename(combined@input.file), ""))
                combined.run@directory = paste("K", k, "/run", re, "/", sep="") 
                dir.create(paste(combined@projDir, combined@lfmmDir, combined.run@directory, 
                    sep = ""), showWarnings = FALSE, recursive = TRUE)
                # lfmmClass.file
                if (to.combine.run@all) {
                    combined.run@lfmmClass.file = paste(tmp, "_r", re ,"_a", to.combine.run@d,".",k, 
                        ".lfmmClass", sep = "")
                } else {
                    combined.run@lfmmClass.file = paste(tmp, "_r", re ,"_s", to.combine.run@d,".",k, 
                        ".lfmmClass", sep = "")
                }
                # run
                combined.run@run = re
                # zscore
                if (to.combine.run@all)
                    combined.run@zscore.file = paste(tmp, "_r", re , "_a", to.combine.run@d, ".",k, ".zscore", sep="")
                else     
                    combined.run@zscore.file = paste(tmp, "_r", re , "_s", to.combine.run@d, ".",k, ".zscore", sep="")
                file.copy(paste(to.combine@projDir, to.combine@lfmmDir, to.combine.run@directory, 
                    to.combine.run@zscore.file, sep = ""), 
                    paste(combined@projDir, combined@lfmmDir, combined.run@directory, 
                    combined.run@zscore.file, sep = "")) 
                # save run
                save.lfmmClass(combined.run, paste(combined@projDir, combined@lfmmDir,
                    combined.run@directory, combined.run@lfmmClass.file, sep = ""))
                # save project
                combined = addRun.lfmmProject(combined, combined.run)
            }
        }
        save.lfmmProject(combined)
    }
)
