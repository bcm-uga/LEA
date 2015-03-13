\name{read.zscore}
\alias{read.zscore}

\title{Read the output files of \code{\link{lfmm}}}
\description{
Read the output file from \code{\link{lfmm}}. This is an internal function.
Zscores of a run can be accessed using the function \code{\link{z.scores}}.
}
\usage{
read.zscore(input.file)
}

\arguments{
\item{input.file}{
a character string containing a path to the output of \code{\link{lfmm}}.
}
}
\value{
\item{R}{
    A matrix containing the \code{\link{lfmm}} results with one line per SNP. 
    The first column is the zscore. The second column is the -log10(p-value).
    The third column is the p-value.
}
}

\author{
    Eric Frichot
}
\seealso{
\code{\link{zscore.format}}
\code{\link{lfmm}}
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
project = lfmm("genotypes.lfmm", "gradients.env", K = 3,
        iterations = 6000, burnin = 3000, project = "new")

res = read.zscore("./genotypes_gradients.lfmm/K3/run1/genotypes_r1_s1.3.zscore")
}


\keyword{read/write}
