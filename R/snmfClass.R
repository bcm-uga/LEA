# create cl sNMF
setClass("snmfClass",
    slots = c(directory="character", 
        snmfClass.file ="character",
        K="integer", run = "integer", 
        CPU = "integer", seed="numeric", alpha = "numeric",
        percentage = "numeric", 
        I = "integer", iterations = "integer",  
        entropy = "logical", tolerance = "numeric",
        crossEntropy = "numeric", ploidy = "integer",
        Q.input.file="character", Q.output.file = "character",
        G.output.file="character")
)

# listMethods

#listMethods_snmfClass <- function()
#{
#c(    "Qvalues",
#    "Gvalues"
#);
#}

# listSlots

#listSlots_snmfClass <- function()
#{
#c(    "directory", "K", "CPU", "seed", "alpha", "missing.data",
#    "snmfClass.file", "percentage ", "I", "iterations", 
#    "entropy", "error", "crossEntropy", "ploidy", "Q.input.file", 
#   "Q.output.file", 
#    "G.output.file")
#}

# .DollarNames.snmfClass <- function(x, pattern) c(listSlots_snmfClass(),
# listMethods_snmfClass())

# $

#setMethod("$", "snmfClass",
#           function(x, name) {
#             if (!(name %in% listMethods_snmfClass() || 
#                name %in% listSlots_snmfClass())) {
#               stop("no $ method for object without attributes")
#         } else if (name %in% listMethods_snmfClass()) {
#        do.call(name, list(x));
#         } else {
#        slot(x, name)
#         }
#           }
#)

# Qvalues

setGeneric("Qvalues", function(object, directory) matrix);
setMethod("Qvalues", "snmfClass",
    function(object, directory) {
        R = as.matrix(read.table(paste(directory, object@directory, 
            object@Q.output.file, sep = "")));
    }
)

# Gvalues

setGeneric("Gvalues", function(object, directory) matrix);
setMethod("Gvalues", "snmfClass",
    function(object, directory) {
        R = as.matrix(read.table(paste(directory, object@directory, 
            object@G.output.file, sep = "")));
    }
)

setGeneric("getCrossEntropy", function(object="snmfClass") vector);
setMethod("getCrossEntropy", "snmfClass",
    function(object) {
        if (object@entropy)
            object@crossEntropy
        else 
            NULL
    }
)
# plot

# display lambda for a value of d, and a Manhattan plot for a value of K. 
#setMethod("plot", "snmfClass",
#          function(x, y, ...){
#        # todo colors
#        barplot(t(x$Qvalues), main="Admixture coefficient plot", ...)
#          }
#)

# show

setMethod("show", "snmfClass",
    function(object) {
        cat("snmf class\n\n")
        cat("file directory:                  ", object@directory, "\n")
        cat("Q output file:                   ", object@Q.output.file, "\n")
        cat("G output file:                   ", object@G.output.file, "\n")
        cat("snmfClass file:                  ", object@snmfClass.file, "\n")
        cat("number of ancestral populations: ", object@K, "\n")
        cat("run number:                      ", object@run, "\n")
        cat("regularization parameter:        ", object@alpha, "\n")
        cat("number of CPUs:                  ", object@CPU, "\n")
        cat("seed:                            ", object@seed, "\n")
        cat("maximal number of iterations:    ", object@iterations, "\n")
        cat("tolerance error:                 ", object@tolerance, "\n")
        cat("Q input file:                    ", object@Q.input.file, "\n")
        if (object@entropy)
            cat("cross-Entropy:                   ", object@crossEntropy,"\n")
        else
            cat("cross-Entropy:                   ", object@entropy,"\n")
    }
)

# summary

#setGeneric("summary", function(object) NULL)
#setMethod("summary", "snmfClass",
#    function(object) {
#        show(object)
#    }
#)

# load

setGeneric("load.snmfClass", function(file="character") attributes("snmfClass"))
setMethod("load.snmfClass", "character",
    function(file) {
        return(dget(file));
    }
)

# save

setGeneric("save.snmfClass", function(x="snmfClass", file="character")NULL)
setMethod("save.snmfClass", signature(x="snmfClass", file="character"),
    function(x, file) {
        dput(x, file) 
    }
)

# remove

setGeneric("remove.snmfClass", function(dir="character", file="character")NULL)
setMethod("remove.snmfClass", signature(dir="character", file="character"),
    function(dir, file) {
        cl = load.snmfClass(paste(dir, file, sep=""))
        file.remove(paste(dir, cl@directory, cl@G.output.file, sep=""))        
        file.remove(paste(dir, cl@directory, cl@Q.output.file, sep=""))        
        file.remove(paste(dir, file, sep=""))
    }
)


