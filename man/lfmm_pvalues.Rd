\name{lfmm.pvalues}
\alias{lfmm.pvalues}

\title{P-values from lfmm runs}
\description{
Returns a vector of p-values computed from a combination of \code{lfmm} runs. For an example, see
\code{\link{lfmm}}.
}
\usage{
lfmm.pvalues (object, genomic.control, lambda, K, d, all, run)
}

\arguments{
\item{object}{
An lfmmProject object.
}
\item{genomic.control}{
A Boolean value. If \code{TRUE}, the p-values are automatically calibrated using genomic control. If \code{FALSE}, the p-values are calculated
by rescaling the chi-squared test statistics using the \code{lambda} parameter.
}
\item{lambda}{A numeric value. The lambda value is used as inflation factor to rescale the chi-squared statistics in the computation of p-values. This option requires that \code{genomic.control} = \code{FALSE}. The default value of \code{lambda} is equal to 1.0 (no rescaling).
}
\item{K}{An integer value. The number of latent factors used in the model.
}
\item{d}{An integer value. Computes the p-values for the d-th covariable in the model.
}
\item{all}{A Boolean value. Each variable is considered separately (Obscure parameter).
}
\item{run}{An integer vector representing a list of runs to be combined in the computation of p-values (by default, all runs).
}
}
\value{
\item{pvalues}{A vector of combined p-values for each locus.}
\item{GIF}{The inflation factor value used for correcting the test statistics.}
}
\author{
    Eric Frichot
    Olivier Francois
}
\seealso{
\code{\link{lfmm.data}}
\code{\link{lfmm}}
}
\examples{
### Example of analyses using lfmm ###

data("tutorial")
# creation of a genotype file, "genotypes.lfmm".
# The data contain 400 SNPs for 50 individuals.
write.lfmm(tutorial.R, "genotypes.lfmm")
# creation of an environmental variable file, "gradient.env".
# The data contain one environmental variable measured for 50 individuals.
write.env(tutorial.C, "gradients.env")

################
# lfmm runs    #
################

# main options, K: (the number of latent factors), 
#           CPU: the number of CPUs.

# runs with K = 3 and 2 repetitions.
# around 15 seconds per run.
project = NULL
project = lfmm("genotypes.lfmm", "gradients.env", K = 3, repetitions = 2,
    iterations = 6000, burnin = 3000, project = "new")

# get adjusted p-values using the genomic control method
p = lfmm.pvalues(project, K = 3)

hist(p$pvalues, col = "yellow3")

# get adjusted p-values using lambda = 0.6
p = lfmm.pvalues(project, genomic.control = FALSE, 
    lambda = 0.6, K = 3)

hist(p$pvalues, col = "yellow3")

}


\keyword{lfmm}
