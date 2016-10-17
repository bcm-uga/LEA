\name{z.scores}
\alias{z.scores}

\title{z-scores from a lfmm run}
\description{
Return the \code{lfmm} output matrix of zscores for the chosen runs with K
latent factors, the d-th variable and the all option. For an example, see
\code{\link{lfmm}}.
}
\usage{
z.scores (object, K, d, all, run)
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
A matrix containing a vector of z-scores for the chosen runs per column. 
}
}
\author{
    Eric Frichot
}
\seealso{
\code{\link{lfmm}}
\code{\link{lfmm.data}}
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

# get the z-scores for all runs for K = 3
z = z.scores(project, K = 3)

# get the z-scores for the 2nd run for K =3
z = z.scores(project, K = 3, run = 2)

# remove
remove.lfmmProject("genotypes_gradients.lfmmProject")
}


\keyword{lfmm}
