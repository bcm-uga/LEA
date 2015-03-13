
# create cl PCA
setClass("pcaProject",
    slots = c(directory="character", 
        n="integer", L="integer", K="integer", 
        center="logical", scale="logical", 
        pcaProject.file ="character",
        input.file ="character",
        eigenvalue.file="character",
        eigenvector.file="character",
        sdev.file="character",
        projection.file="character")
)

# listMethods

listMethods_pcaProject <- function()
{
c(  "eigenvalues",
    "eigenvectors",
    "sdev",
    "projections"
);
}

# listSlots

#listSlots_pcaProject <- function()
#{
#c(    "directory", "n", "L", "K", "center",
#    "scale", "pcaProject_file", "input_file",
#    "eigenvalue_file", "eigenvector_file",
#    "sdev_file", "projection_file")
#}

.DollarNames.pcaProject <- function(x, pattern = "") {
    names <- listMethods_pcaProject()
    grep(pattern, names, value = TRUE )
}


#}listMethods_pcaProject() #c(listSlots_pcaProject(),listMethods_pcaProject())

# [
#setMethod("[", signature(x = "testClass", i = "ANY", j="ANY"),
#         function (x, i, j, ..., drop){
#             print("void function")
#         }
#)

# $

setMethod("$", "pcaProject",
    function(x, name) {

        nb = grep(paste("^",name,sep=""), listMethods_pcaProject(), 
            value = TRUE) 
        if (length(nb) == 0) {
            stop(paste("no $ attribute '", name, "' for the object of class ",
                class(x)[1] ,".",sep=""))
        } else if (length(nb) == 1) {
            do.call(nb, list(x));
        } else {
            stop(paste("Several $ attributes for object of class", class(x)[1],
                "start with", name, ":\n", paste(nb,collapse=" "), sep=""))
        }
    }
)

# eigenvalues

setGeneric("eigenvalues", function(object) matrix);
setMethod("eigenvalues", "pcaProject",
    function (object) {
        as.matrix(read.table(object@eigenvalue.file,sep="/"));
    }
)

# eigenvectors

setGeneric("eigenvectors", function(object) matrix);
setMethod("eigenvectors", "pcaProject",
    function (object) {
        as.matrix(read.table(object@eigenvector.file));
    }
)

# sdev

setGeneric("sdev", function(object) matrix);
setMethod("sdev", "pcaProject",
    function (object) {
        as.matrix(read.table(object@sdev.file));
    }
)

# projections

setGeneric("projections", function(object) matrix);
setMethod("projections", "pcaProject",
    function (object) {
        as.matrix(read.table(object@projection.file));
    }
)

# plot

setMethod("plot", "pcaProject",
    function(x, y, ...){
        plot(eigenvalues(x), ...)
    }
)

# show

setMethod("show", "pcaProject",
    function(object) {
        cat("* pca class *\n\n")
        cat("file directory:                 ", object@directory, "\n")
        cat("input file:                     ", object@input.file, "\n")
        cat("eigenvalue file:                ", 
            basename(object@eigenvalue.file), "\n")
        cat("eigenvector file:               ", 
            basename(object@eigenvector.file), "\n")
        cat("standard deviation file:        ", 
            basename(object@sdev.file), "\n")
        cat("projection file:                ", 
            basename(object@projection.file), "\n")
        cat("pcaProject file:                  ", 
            basename(object@pcaProject.file), "\n")
        cat("number of individuals:          ", object@n, "\n")
        cat("number of loci:                 ", object@L, "\n")
        cat("number of principal components: ", object@K, "\n")
        cat("centered:                       ", object@center, "\n")
        cat("scaled:                         ", object@scale, "\n")
    }
)

# summary

setGeneric("summary", function(object) NULL)
setMethod("summary", "pcaProject",
    function(object) {
        cat("Importance of components:\n");
        rownames = c("Standard deviation", "Proportion of Variance", 
            "Cumulative Proportion");
        colnames = paste("PC", 1:object@n, sep="");
        R = matrix(NA, ncol=object@K, nrow=3, 
            dimnames= list(rownames, colnames));
        R[1,] = sdev(object);
        e = eigenvalues(object);
        R[2,] = e/sum(e);
        R[3,] = cumsum(R[2,]);

        R
    }
)

# load

setGeneric("load.pcaProject", function(file="character") 
    attributes("pcaProject"))
setMethod("load.pcaProject", "character",
    function(file) {
        res = dget(file);

        return(res);
    }
)

# save

setGeneric("save.pcaProject", function(x="pcaProject", file="character") NULL)
setMethod("save.pcaProject", signature(x="pcaProject", file="character"),
    function(x, file) {
        dput(x, file) 
    }
)

# remove

setGeneric("remove.pcaProject", function(file="character")NULL)
setMethod("remove.pcaProject", signature(file="character"),
    function(file) {
        cl = load.pcaProject(file)
        # remove eigenvalues
        file.remove(cl@eigenvalue.file)
        # remove eigenvectors
        file.remove(cl@eigenvector.file)
        # remove standard deviations
        file.remove(cl@sdev.file)
        # remove projections
        file.remove(cl@projection.file)
        # remove tracyWidom file if it exists
        tracywidom.file = setExtension(cl@eigenvalue.file, "tracywidom")
        if (file.exists(tracywidom.file)) {
            file.remove(tracywidom.file)
        }
        # remove directory
        unlink(cl@directory, recursive = TRUE)
        # remove pcaProject file
        file.remove(file)
    }
)

# tracyWidom

setGeneric("tracy.widom", function(object) matrix)
setMethod("tracy.widom", "pcaProject",
    function(object) {

    eigenvalue.input.file = object@eigenvalue.file;

    # output file 
    tracywidom.output.file = setExtension(eigenvalue.input.file, "tracywidom")

    if (!file.exists(tracywidom.output.file)) {

        print("*******************");
        print(" Tracy-Widom tests ");
        print("*******************");

        # test arguments and init
        # input file
        eigenvalue.input.file = test_character("input.file", 
            eigenvalue.input.file, NULL)
        # check extension 
        test_extension(eigenvalue.input.file, "eigenvalues")
        # output file 
        tracywidom.output.file = setExtension(eigenvalue.input.file,
            ".tracywidom")

        .C("R_tracyWidom",
            as.character(eigenvalue.input.file),
            as.character(tracywidom.output.file)
        );
    }
    read.table(tracywidom.output.file, header = TRUE)
        
    }
)
