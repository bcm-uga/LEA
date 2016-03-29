\name{pca}
\alias{pca}
\alias{eigenvalues}
\alias{eigenvectors}
\alias{projections}
\alias{sdev}
\alias{load.pcaProject}
\alias{eigenvalues,pcaProject-method}
\alias{eigenvectors,pcaProject-method}
\alias{import.pcaProject}
\alias{import.pcaProject,character-method}
\alias{export.pcaProject}
\alias{export.pcaProject,character-method}
\alias{projections,pcaProject-method}
\alias{sdev,pcaProject-method}
\alias{load.pcaProject,character-method}
\alias{tracy.widom,pcaProject-method}
\alias{$,pcaProject-method}
\alias{remove.pcaProject}
\alias{remove.pcaProject,character-method}
\alias{plot,pcaProject-method}
\alias{show,pcaProject-method}
\alias{summary,pcaProject-method}

\title{Principal Component Analysis}
\description{
The function \code{pca} performs a Principal Component Analysis 
of a genotypic matrix using the \code{\link{lfmm}}, \code{\link{geno}}, 
\code{\link{ancestrymap}}, \code{\link{ped}} or \code{\link{vcf}} format. 
The function computes eigenvalue, eigenvector, and standard deviation 
for each principal component and the projection of each individual 
on each component. The function pca returns an object of class 
"pcaProject" containing the output data and the input parameters. 
}
\usage{
pca (input.file, K, center = TRUE, scale = FALSE)
}

\arguments{
\item{input.file}{
A character string containg the path to the genotype input file,
a genotypic matrix in the \code{\link{lfmm}} format.
}
\item{K}{
    An integer corresponding to the number of principal components
    calculated. By default, all principal components are calculated.
}
\item{center}{
    A boolean option. If true, the data matrix is centered (default: TRUE).
}
\item{scale}{
    A boolean option. If true, the data matrix is centered and scaled (default:
    FALSE).
}
}
\value{
\code{pca} returns an object of class \code{pcaProject} containing the 
following components:
\item{eigenvalues}{
    The vector of eigenvalues.
}
\item{eigenvectors}{
    The matrix of eigenvectors (one column for each eigenvector). 
}
\item{sdev}{
    The vector of standard deviations. 
}
\item{projections}{
    The matrix of projections (one column for each projection).
}
The following methods can be applied to the object of class pcaProject returned
by \code{pca}: 
\item{plot}{ 
    Plot the eigenvalues.
}
\item{show}{
    Display information about the analysis.
}
\item{summary}{
    Summarize the analysis.
}
\item{tracy.widom}{
    Perform Tracy-Widom tests on the eigenvalues.
}
\item{load.pcaProject(file.pcaProject)}{
    Load the file containing a pcaProject object and return the pcaProject
    object.
}
\item{remove.pcaProject(file.pcaProject)}{
    Erase a \code{pcaProject} object. Caution: All the files associated with 
    the object will be removed.
}
\item{export.pcaProject(file.pcaProject)}{
    Create a zip file containing the full \code{pcaProject} object. It allows to move
    the project to a new directory or a new computer (using import). If you want
    to overwrite an existing export, use the option \code{force == TRUE}.
}
\item{import.pcaProject(file.pcaProject)}{
    Import and load an \code{pcaProject} object from a zip file (made with the export 
    function) into the chosen directory. If you want to overwrite an existing project,
    use the option \code{force == TRUE}.

}
}
\author{
    Eric Frichot
}
\seealso{
\code{\link{lfmm.data}}
\code{\link{snmf}}
\code{\link{lfmm}}
\code{\link{tutorial}}
}
\examples{
# Creation of the genotype file "genotypes.lfmm"
# with 1000 SNPs for 165 individuals.
data("tutorial")
write.lfmm(tutorial.R,"genotypes.lfmm")

#################
# Perform a PCA #
#################

# run of PCA
# Available options, K (the number of PCs calculated), 
#                    center and scale. 
# Creation of   genotypes.pcaProject - the pcaProject object.
#               a directory genotypes.pca containing:
# Create files: genotypes.eigenvalues - eigenvalues,    
#               genotypes.eigenvectors - eigenvectors,
#               genotypes.sdev - standard deviations,
#               genotypes.projections - projections,
# Create a pcaProject object: pc.
pc = pca("genotypes.lfmm", scale = TRUE)

#######################
# Display Information #
#######################

# Display information about the analysis.
show(pc)

# Summarize the analysis.
summary(pc)

#####################
# Graphical outputs #
#####################

par(mfrow=c(2,2))

# Plot eigenvalues.
plot(pc, lwd=5, col="red",xlab=("PCs"),ylab="eigen")

# PC1-PC2 plot.
plot(pc$projections)
# PC3-PC4 plot.
plot(pc$projections[,3:4])

# Plot standard deviations.
plot(pc$sdev)

#############################
# Perform Tracy-Widom tests #
#############################

# Perfom Tracy-Widom tests on all eigenvalues.
# Create file: genotypes.tracyWidom - tracy-widom test information, 
#          in the directory genotypes.pca/.
tw = tracy.widom(pc)

# Plot the percentage of variance explained by each component.
plot(tw$percentage)

# Display the p-values for the Tracy-Widom tests. 
tw$pvalues

##########################
# Manage an pca project #
##########################

# All the file of pca for a given file are 
# automatically saved into a pca project directory and a file.
# The name of the pcaProject file is the same name as 
# the name of the input file with a .pcaProject extension 
# ("genotypes.pcaProject").
# The name of the pcaProject directory is the same name as
# the name of the input file with a .pca extension ("genotypes.pca/")
# There is only one pca Project for each input file including all the runs.

# An pcaProject can be load in a different session.
project = load.pcaProject("genotypes.pcaProject") 

# An pcaProject can be exported to be imported in another directory
# or in another computer
export.pcaProject("genotypes.pcaProject")
dir.create("test/")
#import
newProject = import.pcaProject("genotypes_pcaProject.zip", "test/")
# remove
remove.pcaProject("test/genotypes.pcaProject")

# An pcaProject can be erased.
# Caution: All the files associated with the project will be removed.
remove.pcaProject("genotypes.pcaProject")
}


\keyword{pca}
\keyword{tutorial}
