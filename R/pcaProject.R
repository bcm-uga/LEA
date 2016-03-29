
# create cl PCA
setClass("pcaProject",
    slots = c(projDir="character",
        pcaDir = "character", 
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
        as.matrix(read.table(paste(object@projDir, object@pcaDir,
            object@eigenvalue.file, sep="")));
    }
)

# eigenvectors

setGeneric("eigenvectors", function(object) matrix);
setMethod("eigenvectors", "pcaProject",
    function (object) {
        as.matrix(read.table(paste(object@projDir, object@pcaDir, 
            object@eigenvector.file, sep = "")));
    }
)

# sdev

setGeneric("sdev", function(object) matrix);
setMethod("sdev", "pcaProject",
    function (object) {
        as.matrix(read.table(paste(object@projDir, object@pcaDir, 
            object@sdev.file, sep = "")));
    }
)

# projections

setGeneric("projections", function(object) matrix);
setMethod("projections", "pcaProject",
    function (object) {
        as.matrix(read.table(paste(object@projDir, object@pcaDir, 
            object@projection.file, sep = "")));
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
        cat("project directory:              ", object@projDir, "\n")
        cat("pca result directory:           ", object@pcaDir, "\n")
        cat("input file:                     ", object@input.file, "\n")
        cat("eigenvalue file:                ", object@eigenvalue.file, "\n")
        cat("eigenvector file:               ", object@eigenvector.file, "\n")
        cat("standard deviation file:        ", object@sdev.file, "\n")
        cat("projection file:                ", object@projection.file, "\n")
        cat("pcaProject file:                  ", object@pcaProject.file, "\n")
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
        res = dget(file)
        # remove directory
        unlink(paste(res@projDir, res@pcaDir, sep = ""), recursive = TRUE)
        # remove pcaProject file
        file.remove(file)
    }
)

# tracyWidom

setGeneric("tracy.widom", function(object) matrix)
setMethod("tracy.widom", "pcaProject",
    function(object) {

    eigenvalue.input.file = paste(object@projDir, object@pcaDir, object@eigenvalue.file, sep = "");

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

# export

setGeneric("export.pcaProject", function(file, force) character)
setMethod("export.pcaProject", "character",
    function(file, force = FALSE) {
        object = load.pcaProject(file)
        # entropy
        force = test_logical("force", force, FALSE)

        pathFile = paste(object@projDir, object@pcaProject.file, sep = "")
        zipFile = paste(setExtension(pathFile, ""), "_pcaProject.zip", sep = "")

        if (force == FALSE && file.exists(zipFile)) {
            stop("An export with the same name already exists.\n",
                 "If you want to overwrite it, please use the ",
                 "option 'force == TRUE'\n\n") 
        } else {
            if (file.exists(zipFile))
                file.remove(zipFile)
            curDir = getwd()
            setwd(object@projDir)
            zip(normalizePath(zipFile), c(object@pcaProject.file,
                paste(object@pcaDir, sep = ""), object@input.file))
            setwd(curDir)
            cat("An export of the pca project hase been created: ", currentDir(zipFile), "\n\n") 
        } 
        
    }
)

# import

setGeneric("import.pcaProject", function(zipFile, directory, force) attributes("pcaProject"))
setMethod("import.pcaProject", "character",
    function(zipFile, directory = getwd(), force = FALSE) {
        # force
        force = test_logical("force", force, FALSE)
        # check that no file exists
        tmp = basename(zipFile)
        file = paste(normalizePath(directory), "/", substr(tmp, 1, 
            nchar(tmp) - 15), ".pcaProject", sep = "")
        if (!force && file.exists(file)) {
            stop("A pca project with same name exists in the directory: \n",
                 directory, ".\n",
                 "If you want to overwrite it, please use the ",
                 "option 'force == TRUE'\n\n") 
        }
        # unzip
        unzip(zipFile, exdir = directory)
        cat("The project has been imported into directory, ", directory, "\n\n") 
        res = dget(file);
        res@projDir = paste(normalizePath(directory), "/", sep = "")
        dput(res, file)
        # load the project
        load.pcaProject(file)
    } 
)

