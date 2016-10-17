# create cl lfmm
setClass("lfmmClass",
    slots = c(directory="character", 
        lfmmClass.file ="character",
        K="integer", d = "integer", Niter = "integer", burn = "integer",
        CPU = "integer", seed="numeric", #inflationFactor = "numeric", 
        missing.data = "logical", all="logical", run="integer", 
        epsilon.noise = "numeric", epsilon.b = "numeric",
        random.init = "logical",
        zscore.file="character", 
        deviance="numeric", DIC="numeric"
    )
)

# listMethods

#listMethods_LFMM <- function()
#{
#c(    "zscores",
#    "pvalues",
#    "mlog10pvalues"
#);
#}

# listSlots

#listSlots_LFMM <- function()
#{
#c(    "directory", "K", "Niter", 
#    "burn", "CPU", "seed", "all", "missing.data",
#    "d", "lfmmClass.file", 
#    "noise.epsilon", "b.epsilon", #"inflationFactor",
#    "zscore.file", "dic.file", "epsilon.noise", "epsilon.b", 
#    "random.init", "deviance", "DIC")
#}

# .DollarNames.lfmmClass <- function(x, pattern) 
# c(listSlots_LFMM(),listMethods_LFMM())

# $ 

#setMethod("$", "lfmmClass",
#    function(x, name) {
#        if (!(name %in% listMethods_LFMM() || name %in% listSlots_LFMM())) {
#            stop("no $ method for object without attributes")
#         } else if (name %in% listMethods_LFMM()) {
#        do.call(name, list(x));
#         } else {
#        slot(x, name)
#         }
#           }
#)

# pvalues

setGeneric("pvalues", function(object, directory) matrix);
setMethod("pvalues", "lfmmClass",
    function(object, directory) {
        R = read.zscore(paste(directory, object@directory, object@zscore.file, 
            sep = ""));    
        res = R$pvalues;
    }
)

# mlog10pvalues

setGeneric("mlog10pvalues", function(object, directory) matrix);
setMethod("mlog10pvalues", "lfmmClass",
    function(object, directory) {
        R = read.zscore(paste(directory, object@directory, object@zscore.file, 
            sep = ""));    
        res = R$mlog10pvalues;
    }
)

# zscores

setGeneric("zscores", function(object, directory) matrix);
setMethod("zscores", "lfmmClass",
    function(object, directory) {
        R = read.zscore(paste(directory, object@directory, object@zscore.file, 
            sep = ""));    
        res = R$zscores;
    }
)

# ecdf

#setGeneric("ecdf",  function(object) NULL)
#setMethod("ecdf", "lfmmClass", 
#        function(object){
#        main="Empirical Cumulative Distribution Function"
#        plot.ecdf(pvalues(object), main=main)
#    }
#)

# plot

#setMethod("plot", "lfmmClass",
#          function(x, y, ...){
#            plot(mlog10pvalues(x), main="Manhattan plot", ...)
#          }
#)

# inflationFactorEstimation

#setGeneric("inflationFactorEstimation", function(object) numeric)
#setMethod("inflationFactorEstimation", "lfmmClass",
#    function (object) {
#        pp =  seq(.5,.9, by=0.01)
#        #zs.0 = zscores(object)
#        #p.0 = 1 - pchisq( zs.0^2  , df = 1)
#        print("inflation factor TODO")
#        inflationFactor = 1; #median(quantile(zs.0^2, prob = pp)/qchisq(pp, 
#        df = 1))
#
#        return(inflationFactor)
#    }
#)

# cFDR
#setGeneric("cFDR.lfmmClass",function(object, percentage) vector)
#setMethod("cFDR.lfmmClass", signature=signature(object="lfmmClass", 
#    percentage="numeric"),
#    function(object, percentage) {
#        zs.0 = zscores(object)
#        #inflationFactor = object@inflationFactor
#        p = 1 - pchisq( zs.0^2/inflationFactor , df = 1)
#        k = sum(sort(p) < percentage * (1:length(p)) / length(p))
#        res = which(p <  percentage * k / length(p))
#    } 
#)

# show

setMethod("show", "lfmmClass",
    function(object) {
        cat("* lfmm class *\n\n")
        cat("file directory:                 ", object@directory, "\n")
        cat("lfmmClass file:                 ", object@lfmmClass.file, "\n")
        cat("zscore file:                    ", object@zscore.file,"\n")
        cat("number of latent factors:       ", object@K, "\n")
        cat("run number:                     ", object@run, "\n")
        cat("number of total iterations:     ", object@Niter, "\n")
        cat("number of burnin iterations:    ", object@burn, "\n")
        cat("number of CPUs:                 ", object@CPU, "\n")
        cat("seed:                           ", object@seed, "\n")
        cat("missing data:                   ", object@missing.data, "\n")
        cat("noise epsilon:                  ", object@epsilon.noise, "\n")
        cat("correlation epsilon:            ", object@epsilon.b, "\n")
        cat("random init:                    ", object@random.init, "\n")
        cat("all variables at the same time: ", object@all, "\n")
        if (object@d) 
            cat("Run for variable:               ", object@d, "\n")
        #cat("DIC file:                       ", object@dic.file, "\n")
    }
)

# summary

#setGeneric("summary", function(object) NULL)
#setMethod("summary", "lfmmClass",
#    function(object) {
#        print("TODO")
#        percentage = c(.9,.8,.7,.6,.5,.4,.3,.2,.1,.05,.02,.01)
#        colnames = paste(percentage *100, "%")    
#        rownames = "proportion of correlations"
#        res = matrix(NA,ncol=length(percentage), nrow=1, 
#            dimnames= list(rownames, colnames)) 
#        for (p in 1:length(percentage))
#            res[p] = length(cFDR.lfmmClass(object, percentage[p]))/object@L
#        return(list(cFDR=res));
#    }
#)

# load

setGeneric("load.lfmmClass", function(file="character") attributes("lfmmClass"))
setMethod("load.lfmmClass", "character",
    function(file) {
        return(dget(file));
    }
)

# save

setGeneric("save.lfmmClass", function(x="lfmmClass", file="character") NULL)
setMethod("save.lfmmClass", signature(x="lfmmClass", file="character"),
    function(x, file) {
        dput(x, file) 
    }
)

# remove

setGeneric("remove.lfmmClass", function(dir="character", file="character")NULL)
setMethod("remove.lfmmClass", signature(dir="character", file="character"),
    function(dir, file) {
        cl = load.lfmmClass(paste(dir, file, sep = ""))
        file.remove(paste(dir, cl@directory, cl@zscore.file, sep = ""))
        file.remove(file)
    }
)

