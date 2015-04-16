\name{mlog10p.values}
\alias{mlog10p.values}

\title{-log10(p-values) from a lfmm run}
\description{
Return the \code{lfmm} output matrix of -log10(p-values) for the chosen runs 
with K latent factors, the d-th variable and the all option. For an example, 
see \code{\link{lfmm}}.
}
\usage{
mlog10p.values (object, K, d, all, run)
}

\arguments{
\item{object}{
A lfmmProject object.
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
A matrix containing a vector of -log10(p-values) for the chosen runs per 
column. 
}
}
\author{
    Eric Frichot
}
\seealso{
\code{\link{lfmm.data}}
\code{\link{lfmm}}
\code{\link{p.values}}
\code{\link{z.scores}}
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

# get the -log10(p-values) for all runs for K = 3
mp = mlog10p.values(project, K = 3)

# get the -log10(p-values) for the 2nd run for K =3
mp = mlog10p.values(project, K = 3, run = 2)
}


\keyword{lfmm}
