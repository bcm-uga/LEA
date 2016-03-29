\name{adjusted.pvalues}
\alias{adjusted.pvalues}

\title{adjusted p-values from a lfmm run}
\description{
Return the \code{lfmm} output vector of adjusted p-values and the genomic
inflation factor using the genomic control method or the lambda inflation 
factor parameter for the chosen runs with K fatent factors, the d-th 
variable and the all option. For an example, see
\code{\link{lfmm}}.
}
\usage{
adjusted.pvalues (object, genomic.control, lambda, K, d, all, run)
}

\arguments{
\item{object}{
A lfmmProject object.
}
\item{genomic.control}{
A boolean option. If true, the p-values are automatically calibrated
using the genomic control method. If false, the p-values are calculated
using the \code{lambda} inflation factor parameter.
}
\item{lambda}{
the lambda inflation factor used to calibrate the p-value if 
\code{genomic.control} = FALSE (default: 1.0).
}
\item{K}{
The number of latent factors.
}
\item{d}{
The d-th variable.
}
\item{all}{
A Boolean option. If true, the run with all variables at the same time. If 
false, the runs with each variable separately.
}
\item{run}{
A list of chosen runs.
}
}
\value{
\item{res}{
A matrix containing a vector of p.values for the chosen runs per column. 
}
}
\author{
    Eric Frichot
}
\seealso{
\code{\link{lfmm.data}}
\code{\link{lfmm}}
\code{\link{p.values}}
\code{\link{mlog10p.values}}
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

# Toy runs with K = 3 and 2 repetitions.
# around 15 seconds per run.
project = NULL
project = lfmm("genotypes.lfmm", "gradients.env", K = 3, repetitions = 2,
    iterations = 6000, burnin = 3000, project = "new")

# get the adjusted p-values using the genomic control method
res = adjusted.pvalues(project, K = 3)

hist(res$p.values, col = "yellow3")

# get the adjusted p-values with the genomic inflatino factor
res = adjusted.pvalues(project, genomic.control = FALSE, 
    lambda = res$genomic.inflation.factor, K = 3)

hist(res$p.values, col = "yellow3")

}


\keyword{lfmm}
