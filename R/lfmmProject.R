# create Class LFMM Project
setClass("lfmmProject",
    slots = c(lfmmProject.file = "character", directory="character", 
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
            run@lfmmClass.file);

        return(project)
    }
)


# listMethods

#listMethods_lfmmProject <- function()
#{
#}

# listSlots

#listSlots_lfmmProject <- function()
#{
#c(    "snmfProject_file", "input_file", "runs", "K", "d", "all", 
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
        cat("directory:                            ", object@directory, "\n")
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
            res[,i] = mlog10pvalues(object@runs[[r[run[i]]]])        
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
            res[,i] = pvalues(object@runs[[r[run[i]]]])        
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
            res[,i] = zscores(object@runs[[r[run[i]]]])        
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
                res@runs[[r]] = load.lfmmClass(res@lfmmClass.files[r])
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
                save.lfmmClass(object@runs[[r]], object@lfmmClass.files[r])
            }
        }
        object@runs = list()
        dput(object, file) 
        cat("The project is saved into :\n",currentDir(file),"\n\n");
        cat("To load the project, use:\n project = load.lfmmProject(\"",
            currentDir(file),"\")\n\n",sep="");
        cat("To remove the project, use:\n remove.lfmmProject(\"",
            currentDir(file), "\")\n\n", sep="")

        object@lfmmProject.file;
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
        unlink(res@directory, recursive = TRUE)
        file.remove(file)
    }
)

#setGeneric("inflationFactorEstimation", function(object, ...) vector)
#setMethod("inflationFactorEstimation", "lfmmProject",
#        function(object, K, run, d, all , type) {
#    }
#)
