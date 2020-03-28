\name{lfmm}
\alias{lfmm}
\alias{load.lfmmProject}
\alias{remove.lfmmProject}
\alias{load.lfmmProject,character-method}
\alias{remove.lfmmProject,character-method}
\alias{import.lfmmProject}
\alias{import.lfmmProject,character-method}
\alias{export.lfmmProject}
\alias{export.lfmmProject,character-method}
\alias{combine.lfmmProject}
\alias{combine.lfmmProject,character,character-method}
\alias{show,lfmmProject-method}
\alias{show,lfmmClass-method}
\alias{summary,lfmmProject-method}
\alias{z.scores,lfmmProject-method}
\alias{p.values,lfmmProject-method}
\alias{adjusted.pvalues,lfmmProject-method}
\alias{mlog10p.values,lfmmProject-method}
\alias{plot,lfmmProject-method}
\alias{lfmm.pvalues,lfmmProject-method}


\title{Fitting Latent Factor Mixed Models (MCMC algorithm)}
\description{Latent Factor Mixed Models (LFMMs) are factor regression models in which the response variable is a genotypic matrix, and the explanatory variables are environmental measures of ecological interest or trait values. The \code{\link{lfmm}} function estimates latent factors and effect sizes based on an MCMC algorithm. The resulting object can be used in the function \code{\link{lfmm.pvalues}} to identify genetic polymorphisms exhibiting association with ecological gradients or phenotypes, while correcting for unobserved confounders. An exact and computationally efficient least-squares method is implemented in the function \code{\link{lfmm2}}.
}

\usage{
lfmm(input.file, environment.file, K, 
    project = "continue", 
    d = 0, all = FALSE, 
    missing.data = FALSE, CPU = 1, 
    iterations = 10000, burnin = 5000, 
    seed = -1, repetitions = 1, 
    epsilon.noise = 1e-3, epsilon.b = 1000, 
    random.init = TRUE)

}

\arguments{
\item{input.file}{
    A character string containing a path to the input file,
    a genotypic matrix in the \code{\link{lfmm}{lfmm_format}} format. The matrix must not contain missing values. See \code{\link{impute}} for completion based on nonnegative matrix factorization.
}
\item{environment.file}{
    A character string containing a path to the environmental file,
    an environmental data matrix in the \code{\link{env}} format.
}
\item{K}{
    An integer corresponding to the number of latent factors. 
}
\item{project}{
    A character string among "continue", "new", and "force". If "continue",
    the results are stored in the current project. If "new", the current project
    is removed and a new project is created. If "force", the 
    results are stored in the current project even if the input file has been 
    modified since the creation of the project.
}
\item{d}{
    An integer corresponding to the fit of an \code{lfmm} model with the 
    \code{d}-th variable only from \code{environment.file}. By default (if NULL
    and all are FALSE), \code{lfmm} fits each variable from 
    \code{environment.file} sequentially and independently.
}
\item{all}{
    A Boolean option. If TRUE, \code{lfmm} fits all variables from the
    \code{environment.file} at the same time. This option is not compatible 
    with the \code{d} option. 
}
\item{missing.data}{
    A Boolean option. If TRUE, the \code{input.file} contains missing genotypes. Caution: \code{lfmm} requires imputed genotype matrices. See 
    \code{\link{impute}}.
}
#ifndef windows
\item{CPU}{
    A number of CPUs to run the parallel version of the algorithm.
    By default, the number of CPUs is 1.
}
#endif
#ifdef windows
\item{CPU}{
    The windows version is not parallelized. By consequence, CPU is always set to 1.
}
#endif
\item{iterations}{
    The total number of cycles for the Gibbs Sampling algorithm.
}
\item{burnin}{
    The burnin number of cycles for the Gibbs Sampling algorithm.
}
\item{seed}{
    A seed to initialize the random number generator. By default, the seed is 
    randomly chosen. The seed is initialized in each run of the program. 
    For modifying the default setting, provide one seed per run.
}
\item{repetitions}{
    A number of replicate runs for the Gibbs Sampler algorithm.
}
\item{epsilon.noise}{
    A prior parameter for variances.
}
\item{epsilon.b}{
    A prior parameter for the variance of correlation coefficients.
}
\item{random.init}{
    A Boolean option. If TRUE, the Gibbs Sampler is initiliazed randomly. 
    Otherwise, it is initialized with zero values.
}
}
\value{
\code{lfmm} returns an object of class \code{lfmmProject}.

The following methods can be applied to an object of class \code{lfmmProject}:

\item{show}{
    Display information about all analyses.
}
\item{summary}{
    Summarize analyses.
}
\item{\code{\link{z.scores}}}{
    Return the \code{lfmm} output vector of z.scores for some runs.
}
\item{\code{\link{lfmm.pvalues}}}{
    Return the vector of adjusted p-values for a combination of runs with K
    latent factors, and for the d-th predictor.
}
\item{load.lfmmProject (file = "character") }{
    Load the file containing an lfmmProject objet and show the
    object.
}
\item{remove.lfmmProject (file = "character")}{
    Erase a \code{lfmmProject} object. Caution: All the files associated with 
    the object will be removed.
}
\item{export.lfmmProject(file.lfmmProject)}{
    Create a zip file containing the full \code{lfmmProject} object. It allows users to move
    the project to a new directory or a new computer (using import). If you want
    to overwrite an existing export, use the option \code{force == TRUE}.
}
\item{import.lfmmProject(file.lfmmProject)}{
    Import and load an \code{lfmmProject} object from a zip file (made with the export 
    function) into the chosen directory. If you want to overwrite an existing project,
    use the option \code{force == TRUE}.
}
\item{combine.lfmmProject(file.lfmmProject, toCombine.lfmmProject)}{
    Combine \code{to.Combine.lfmmProject} into \code{file.lfmmProject}.
    Caution: Only projects with runs coming from the same input file can be combined.
    If the same input file has different names in the two projects, use the option
    \code{force == TRUE}.
}
}
\author{
    Eric Frichot
    Olivier Francois
}
\seealso{
\code{\link{lfmm.data}}
\code{\link{z.scores}}
\code{\link{lfmm.pvalues}}
\code{\link{pca}}
\code{\link{lfmm}}
\link{tutorial}
}
\references{
Frichot E, Schoville SD, Bouchard G, Francois O. (2013). \emph{Testing for 
associations between loci and environmental gradients using latent factor mixed
models}. Molecular biology and evolution, 30(7), 1687-1699.
}
\examples{
### Example of analysis using lfmm ###

data("tutorial")
# creation of a genotype file: genotypes.lfmm.
# The file contains 400 SNPs for 50 individuals.
write.lfmm(tutorial.R, "genotypes.lfmm")

# Creation of a phenotype/environment file: gradient.env.
# One environmental predictor for 40 individuals.
write.env(tutorial.C, "gradients.env")

################
# running lfmm #
################

# main options, K: (the number of latent factors), 
#           CPU: the number of CPUs.

# Runs with K = 6 and 5 repetitions.
# runs with 6000 iterations 
# including 3000 iterations for burnin.
# Around 30 seconds per run.
project = lfmm( "genotypes.lfmm", 
                "gradients.env", 
                 K = 6, 
                 repetitions = 5, 
                 project = "new")

# get adjusted p-values using all runs
pv = lfmm.pvalues(project, K = 6)

# Evaluate FDR and POWER (TPR)
for (alpha in c(.05,.1,.15,.2)) {
    # expected FDR
    print(paste("expected FDR:", alpha))
    L = length(pv$pvalues)
    # Benjamini-Hochberg's mehod for an expected FDR = alpha.
    w = which(sort(pv$pvalues) < alpha * (1:L)/L)
    candidates = order(pv$pvalues)[w]

    # estimated FDR and True Positive Rate
    # The targets SNPs are loci 351 to 400
    Lc = length(candidates)
    estimated.FDR = length(which(candidates <= 350))/Lc
    estimated.TPR = length(which(candidates > 350))/50
    print(paste("FDR:", estimated.FDR, "True Positive Rate:", estimated.TPR))
}

###################
# Post-treatments #
###################

# show the project
show(project)

# summary of the project
summary(project)

# get the z-scores for the 2nd run for K = 6
z = z.scores(project, K = 6, run = 2)

# get the p-values for K = 6 and run 2
p = lfmm.pvalues(project, K = 6, run = 2)

##########################
# Manage an lfmm project #
##########################

# All the runs of lfmm for a given file are 
# automatically saved into an lfmm project directory and a file.
# The name of the lfmmProject file is the concatenation of 
# the name of the input file and the environment file 
# with a .lfmmProject extension ("genotypes_gradient.lfmmProject").
# The name of the lfmmProject directory is the same name as
# the lfmmProject file with a .lfmm extension ("genotypes_gradient.lfmm/")
# There is a unique lfmm Project for each input file.

# An lfmmProject can be loaded in an R session as follows
project = load.lfmmProject("genotypes_gradients.lfmmProject")

# An lfmmProject can be exported to be imported in another directory
# or in another computer as follows
export.lfmmProject("genotypes_gradients.lfmmProject")

#ifndef windows
dir.create("test", showWarnings = TRUE)
#import
newProject = import.lfmmProject("genotypes_gradients_lfmmProject.zip", "test")

# combine projects
combinedProject <- combine.lfmmProject(
                  "genotypes_gradients.lfmmProject", 
                  "test/genotypes_gradients.lfmmProject"
                  )

# remove
remove.lfmmProject("test/genotypes_gradients.lfmmProject")
#endif

#ifdef windows
# remove
remove.lfmmProject("genotypes_gradients.lfmmProject")

#import
newProject = import.lfmmProject("genotypes_gradients_lfmmProject.zip")
#endif


# An lfmmProject can be removed as follows.
# Caution: All the files associated with the project will be removed.
remove.lfmmProject("genotypes_gradients.lfmmProject")

}


\keyword{lfmm}
\keyword{lfmm2}
\keyword{tutorial}
