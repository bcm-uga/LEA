#create Class sNMF Project
setClass("snmfProject",
    slots = c(snmfProject.file = "character", directory = "character", 
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
        project@snmfClass.files = c(project@snmfClass.files, run@snmfClass.file)

        return(project)
    }
)

# getRuns

setGeneric("getRuns.snmfProject", function(object, ...) 
    standardGeneric("getRuns.snmfProject"));
setMethod("getRuns.snmfProject", "snmfProject",
    function(object, K) {
        # check of K
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

        plot(s$crossEntropy[1,], ylab="Minimal Cross-Entropy", 
            xlab="Number of ancestral populations", ...)
    }
)

# G

setGeneric("G",  function(object, K, run) matrix)
setMethod("G", "snmfProject",
    function(object, K, run) {
        
        # check of K
        if (missing(K)) {
            # if only one, that is the one
            if (length(unique(object@K)) == 1) {
                K = unique(object@K)
            } else {
                stop("Please, choose a value of K among: ", 
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
            
        R = Gvalues(object@runs[[r[run]]])

        return(R)
    }
)

# Q

setGeneric("Q",  function(object, K, run) matrix)
setMethod("Q", "snmfProject",
    function(object, K, run) {
        
        # check of K
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
            
        R = Qvalues(object@runs[[r[run]]])

        return(R)
    }
)

# crossEntropy

setGeneric("cross.entropy", function(object, K, run) vector)
setMethod("cross.entropy", "snmfProject",
    function(object, K, run) {

        # check of K
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

# show

setMethod("show", "snmfProject",
    function(object) {
        cat("snmf Project\n\n")
        cat("snmfProject file:                ", object@snmfProject.file, "\n")
        cat("directory:                       ", object@directory, "\n")
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
                res@runs[[r]] = load.snmfClass(res@snmfClass.files[r])
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
                save.snmfClass(object@runs[[r]], object@snmfClass.files[r])
            }
        }
        object@runs = list()
        dput(object, file) 
        cat("The project is saved into :\n",currentDir(file),"\n\n");
        cat("To load the project, use:\n project = load.snmfProject(\"",
            currentDir(file),"\")\n\n",sep="");
        cat("To remove the project, use:\n remove.snmfProject(\"",
            currentDir(file),"\")\n\n",sep="");
   
        object@snmfProject.file;
    }
)

# remove

setGeneric("remove.snmfProject", function(file="character") NULL)
setMethod("remove.snmfProject", "character",
    function(file) {
        res = dget(file);
        unlink(res@directory, recursive = TRUE)
        file.remove(file)
    }
)
