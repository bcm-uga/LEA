\name{snmf}
\alias{snmf}
\alias{plot,snmfProject-method}
\alias{cross.entropy,snmfProject-method}
\alias{load.snmfProject}
\alias{load.snmfProject,character-method}
\alias{remove.snmfProject}
\alias{remove.snmfProject,character-method}
\alias{import.snmfProject}
\alias{import.snmfProject,character-method}
\alias{export.snmfProject}
\alias{export.snmfProject,character-method}
\alias{combine.snmfProject}
\alias{combine.snmfProject,character,character-method}
\alias{show,snmfProject-method}
\alias{show,snmfClass-method}
\alias{summary,snmfProject-method}
\alias{Q,snmfProject-method}
\alias{G,snmfProject-method}

\title{Estimates individual ancestry coefficients and ancestral allele 
frequencies.}
\description{
\code{\link{snmf}} estimates admixture coefficients using sparse Non-Negative 
    Matrix Factorization algorithms, and provide STRUCTURE-like outputs.
}
\usage{
snmf (input.file, K, 
    project = "continue",
    repetitions = 1, CPU = 1, 
    alpha = 10, tolerance = 0.00001, entropy = FALSE, percentage = 0.05,
    I, iterations = 200, ploidy = 2, seed = -1, Q.input.file)
}

\arguments{
\item{input.file}{
    A character string containing a the path to the input file, a genotypic 
    matrix in the \code{\link{geno}} format.
}
\item{K}{
    An integer vector corresponding to the number of ancestral populations for
    which the snmf algorithm estimates have to be calculated.
}
\item{project}{
    A character string among "continue", "new", and "force". If "continue",
    the results are stored in the current project. If "new", the current 
    project is removed and a new one is created to store the result. If 
    "force", the results are stored in the current project even if the input 
    file has been modified since the creation of the project.
}
\item{repetitions}{
    An integer corresponding with the number of repetitions for each value of 
    \code{K}.
}
#ifndef windows
\item{CPU}{
    A number of CPUs to run the parallel version of the 
    algorithm. By default, the number of CPUs is 1.
}
#endif
#ifdef windows
\item{CPU}{
    The windows version is not parallelized. By consequence, 
    CPU is always 1.
}
#endif
\item{alpha}{
    A numeric value corresponding to the snmf regularization parameter. 
    The results can depend on the value of this parameter, especially for 
    small data sets. 
}
\item{tolerance}{
    A numeric value for the tolerance error.
}
\item{entropy}{
    A boolean value. If true, the cross-entropy criterion is calculated 
    (see \code{\link{create.dataset}} and 
    \code{\link{cross.entropy.estimation}}).
}
\item{percentage}{
    A numeric value between 0 and 1 containing the percentage of 
    masked genotypes when computing the cross-entropy 
    criterion. This option applies only if \code{entropy == TRUE} 
    (see \code{\link{cross.entropy}}).
}
\item{I}{
    The number of SNPs to initialize the algorithm. It starts the algorithm 
    with a run of snmf using a subset of nb.SNPs random SNPs. If this option 
    is set with nb.SNPs, the number of randomly chosen SNPs is the minimum 
    between 10000 and 10 \% of all SNPs. This option can considerably speeds 
    up snmf estimation for very large data sets.
}
\item{iterations}{
    An integer for the maximum number of iterations in algorithm. 
}
\item{ploidy}{
    1 if haploid, 2 if diploid, n if n-ploid.
}
\item{seed}{
    A seed to initialize the random number generator. 
    By default, the seed is randomly chosen.
}
\item{Q.input.file}{
    A character string containing a path to an initialization file for Q, 
    the individual admixture coefficient matrix.
}
}
\value{
\code{snmf} returns an object of class \code{snmfProject}.

The following methods can be applied to the object of class {snmfProject}:
\item{plot}{
    Plot the minimal cross-entropy in function of K.
}
\item{show}{
    Display information about the analyses.
}
\item{summary}{
    Summarize the analyses.
}
\item{\code{\link{Q}}}{
    Return the admixture coefficient matrix for the chosen run with K
    ancestral populations. 
}
\item{\code{\link{G}}}{
    Return the ancestral allele frequency matrix for the chosen run with K
    ancestral populations. 
}
\item{\code{\link{cross.entropy}}}{
    Return the cross-entropy criterion for the chosen runs with K
    ancestral populations. 
}
\item{load.snmfProject(file.snmfProject)}{
    Load the file containing an snmfProject objet and return the snmfProject 
    object.
}
\item{remove.snmfProject(file.snmfProject)}{
    Erase a \code{snmfProject} object. Caution: All the files associated with 
    the object will be removed.
}
\item{export.snmfProject(file.snmfProject)}{
    Create a zip file containing the full \code{snmfProject} object. It allows to move
    the project to a new directory or a new computer (using import). If you want
    to overwrite an existing export, use the option \code{force == TRUE}.
}
\item{import.snmfProject(file.snmfProject)}{
    Import and load an \code{snmfProject} object from a zip file (made with the export 
    function) into the chosen directory. If you want to overwrite an existing project,
    use the option \code{force == TRUE}.
}
\item{combine.snmfProject(file.snmfProject, toCombine.snmfProject)}{
    Combine \code{to.Combine.snmfProject} into \code{file.snmfProject}.
    Caution: Only projects with runs coming from the same input file can be combined.
    If the same input file has different names in the two projects, use the option
    \code{force == TRUE}.
}
}
\author{
    Eric Frichot
}
\seealso{
\code{\link{geno}}
\code{\link{pca}}
\code{\link{lfmm}}
\code{\link{tutorial}}
}
\references{
Frichot E, Mathieu F, Trouillon T, Bouchard G, Francois O. (2014). \emph{Fast 
and Efficient Estimation of Individual Ancestry Coefficients}. Genetics, 
194(4): 973--983.
}
\examples{
### Example of analyses using snmf ###

# creation of the genotype file, genotypes.geno.
# It contains 400 SNPs for 50 individuals.
data("tutorial")
write.geno(tutorial.R, "genotypes.geno")

################
# runs of snmf #
################

# main options, K: (the number of ancestral populations), 
#        entropy: calculate the cross-entropy criterion, 
#        CPU: the number of CPUs.

# Runs with K between 1 and 5 with cross-entropy and 2 repetitions.
project = NULL
project = snmf("genotypes.geno", K=1:10, entropy = TRUE, repetitions = 10,
    project = "new")

# plot cross-entropy criterion of all runs of the project
plot(project, lwd = 5, col = "red", pch=1)

# get the cross-entropy of each run for K = 4
ce = cross.entropy(project, K = 4)

# select the run with the lowest cross-entropy
best = which.min(ce)

# plot the best run for K = 4 (ancestry coefficients).
barplot(t(Q(project, K = 4, run = best)))

###################
# Post-treatments #
###################

# show the project
show(project)

# summary of the project
summary(project)

# get the cross-entropy for all runs for K = 4
ce = cross.entropy(project, K = 4)

# get the cross-entropy for the 2nd run for K = 4
ce = cross.entropy(project, K = 4, run = 2)

# get the ancestral genotype frequency matrix, G, for the 2nd run for K = 4. 
res = G(project, K = 4, run = 2)

#############################
# Advanced snmf run options #
#############################

# Q.input.file: init a run with a given ancestry coefficient matrix Q. 
# Here, it is initialized with the Q matrix from the first run with K=4
project = snmf("genotypes.geno", K = 4, 
    Q.input.file = "./genotypes.snmf/K4/run1/genotypes_r1.4.Q")

# I: init the Q matrix of a run from a smaller run with 100 randomly chosen 
# SNPs. 
project = snmf("genotypes.geno", K = 4, I = 100)

# CPU: run snmf with 2 CPUs. 
project = snmf("genotypes.geno", K = 4, CPU=2)

# percentage: run snmf and calculate the cross-entropy criterion with 10% of 
# masked genotypes, instead of 5% of masked genotypes. 
project = snmf("genotypes.geno", K = 4, entropy= TRUE, percentage = 0.1)

# seed: choose the seed to init the randomization.
project = snmf("genotypes.geno", K = 4, seed=42)

# alpha: choose the regularization parameter.
project = snmf("genotypes.geno", K = 4, alpha = 100)

# tolerance: choose the tolerance parameter.
project = snmf("genotypes.geno", K = 4, tolerance = 0.0001)

##########################
# Manage an snmf project #
##########################

# All the runs of snmf for a given file are 
# automatically saved into a snmf project directory and a file.
# The name of the snmfProject file is the same name as 
# the name of the input file with a .snmfProject extension 
# ("genotypes.snmfProject").
# The name of the snmfProject directory is the same name as
# the name of the input file with a .snmf extension ("genotypes.snmf/")
# There is only one snmf Project for each input file including all the runs.

# An snmfProject can be load in a different session.
project = load.snmfProject("genotypes.snmfProject") 

# An snmfProject can be exported to be imported in another directory
# or in another computer
export.snmfProject("genotypes.snmfProject")
dir.create("test/")
#import
newProject = import.snmfProject("genotypes_snmfProject.zip", "test/")
# combine projects
combinedProject = combine.snmfProject("genotypes.snmfProject", "test/genotypes.snmfProject")
# remove
remove.snmfProject("test/genotypes.snmfProject")

# An snmfProject can be erased.
# Caution: All the files associated with the project will be removed.
remove.snmfProject("genotypes.snmfProject")

}


\keyword{snmf}
\keyword{tutorial}
