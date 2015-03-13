\name{lfmm}
\alias{lfmm}
\alias{load.lfmmProject}
\alias{remove.lfmmProject}
\alias{load.lfmmProject,character-method}
\alias{remove.lfmmProject,character-method}
\alias{show,lfmmProject-method}
\alias{show,lfmmClass-method}
\alias{summary,lfmmProject-method}
\alias{z.scores,lfmmProject-method}
\alias{p.values,lfmmProject-method}
\alias{mlog10p.values,lfmmProject-method}
\alias{plot,lfmmProject-method}


\title{Fitting Latent Factor Mixed Models}
\description{
\code{\link{lfmm}} is used to fit Latent Factor Mixed Models. The goal of 
\code{lfmm} is to identify genetic polymorphisms that exhibit high correlation
with some environmental gradient or with the variables used as proxies for 
ecological pressures.
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
    a genotypic matrix in the \code{\link{lfmm}{lfmm_fomat}} format.
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
    is removed and a new one is created to store the result. If "force", the 
    results are stored in the current project even if the input file has been 
    modified since the creation of the project.
}
\item{d}{
    An integer corresponding to the fit of \code{lfmm} model with the 
    \code{d}-th variable only from \code{environment.file}. By default (if NULL
    and all is FALSE), fit \code{lfmm} with each variable from 
    \code{environment.file} sequentially and independently.
}
\item{all}{
    A boolean option. If true, fit \code{lfmm} with all variables from 
    \code{environment.file} at the same time. This option is not compatible 
    with the \code{d} option. 
}
\item{missing.data}{
    A boolean option. If true, the \code{input.file} contains missing genotypes.
}
#ifndef windows
\item{CPU}{
    A number of CPUs to run the parallel version of the algorithm.
    By default, the number of CPUs is 1.
}
#endif
#ifdef windows
\item{CPU}{
    The windows version is not parallelized. By consequence, CPU is always set
    to 1.
}
#endif
\item{iterations}{
    The total number of iterations in the Gibbs Sampling algorithm.
}
\item{burnin}{
    The burnin number of iterations in the Gibbs Sampling algorithm.
}
\item{seed}{
    A seed to initialize the random number generator. By default, the seed is 
    randomly chosen. The seed is initialized at each repetition. If you want
    to set a seed, please provide a seed per repetition.
}
\item{repetitions}{
    The number of repetitions of each run.
}
\item{epsilon.noise}{
    Prior on the different variances.
}
\item{epsilon.b}{
    Prior on the variance of the correlation coefficients.
}
\item{random.init}{
    A boolean option. If true, the Gibbs Sampler is initiliazed randomly. 
    Otherwise, it is initialized with zeros.
}
}
\value{
\code{lfmm} returns an object of class \code{lfmmProject}.

The following methods can be applied to the object of class \code{lfmmProject}:



\item{show}{
    Display information about the analyses.
}
\item{summary}{
    Summarize the analyses.
}
\item{\code{\link{z.scores}}}{
    Return the \code{lfmm} output vector of zscores for the chosen runs with K
    latent factors, the d-th variable and the all option.
}
\item{\code{\link{p.values}}}{
    Return the \code{lfmm} output vector of p-values for the chosen runs with K
    latent factors, the d-th variable and the all option.
}
\item{\code{\link{mlog10p.values}}}{
    Return the \code{lfmm} output vector of -log10(p-values) for the chosen 
    runs with K latent factors, the d-th variable and the all option.
}
\item{load.lfmmProject (file = "character") }{
    Load the file containing an lfmmProject objet and return the lfmmProject 
    object.
}
\item{remove.lfmmProject (file = "character")}{
    Erase a \code{lfmmProject} object. Caution: All the files associated with 
    the object will be removed.
}
}
\author{
    Eric Frichot
}
\seealso{
\code{\link{lfmm.data}}
\code{\link{z.scores}}
\code{\link{p.values}}
\code{\link{mlog10p.values}}
\code{\link{pca}}
\code{\link{snmf}}
\link{tutorial}
}
\references{
Frichot E, Schoville SD, Bouchard G, Francois O. (2013). \emph{Testing for 
associations between loci and environmental gradients using latent factor mixed
models}. Molecular biology and evolution, 30(7), 1687-1699.
}
\examples{
### Example of analyses using lfmm ###

data("tutorial")
# creation of the genotype file, genotypes.lfmm.
# It contains 400 SNPs for 50 individuals.
write.lfmm(tutorial.R, "genotypes.lfmm")
# creation of the environment file, gradient.env.
# It contains 1 environmental variable for 40 individuals.
write.env(tutorial.C, "gradients.env")

################
# runs of lfmm #
################

# main options, K: (the number of latent factors), 
#           CPU: the number of CPUs.

# Runs with K = 9 and 5 repetitions.
# The runs are composed of 6000 iterations including 3000 iterations
# for burnin.
# around 30 seconds per run.
project = NULL
project = lfmm("genotypes.lfmm", "gradients.env", K = 6, repetitions = 5, 
        project = "new")

# get the zscores of each run for K = 6
zs = z.scores(project, K = 6)

# Combine the z-scores using the Stouffer method
zs.stouffer = apply(zs, MARGIN = 1, median)

# calculate the inflation factor
lambda = median(zs.stouffer^2)/.456

# calculate adjusted p-values
cp.values = pchisq(zs.stouffer^2/lambda, df = 1, lower = FALSE)

for (alpha in c(.05,.1,.15,.2)) {
    # expected FDR
    print(paste("expected FDR:", alpha))
    L = length(cp.values)
    # return a list of candidates with an expected FDR of alpha.
    w = which(sort(cp.values) < alpha * (1:L) / L)
    candidates = order(cp.values)[w]

    # estimated FDR and True Positif
    estimated.FDR = length(which(candidates <= 350))/length(candidates)
    estimated.TP = length(which(candidates > 350))/50
    print(paste("FDR:", estimated.FDR, "True Positive:", estimated.TP))
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

# get the p-values for the 2nd run for K = 6
p = p.values(project, K = 6, run = 2)

# get the -log10(p-values) for the 2nd run for K = 6
mp = mlog10p.values(project, K = 6, run = 2)

##########################
# Manage an lfmm project #
##########################

# All the runs of lfmm for a given file are 
# automatically saved into a lfmm project directory and a file.
# The name of the lfmmProject file is a combination of 
# the name of the input file and the environment file 
# with a .lfmmProject extension ("genotypes_gradient.lfmmProject").
# The name of the lfmmProject directory is the same name as
# the lfmmProject file with a .lfmm extension ("genotypes_gradient.lfmm/")
# There is only one lfmm Project for each input file including all the runs.

# An lfmmProject can be load in a different session.
project = load.lfmmProject("genotypes_gradients.lfmmProject")

# An lfmmProject can be erased.
# Caution: All the files associated with the project will be removed.
remove.lfmmProject("genotypes_gradients.lfmmProject")

}


\keyword{lfmm}
\keyword{tutorial}
