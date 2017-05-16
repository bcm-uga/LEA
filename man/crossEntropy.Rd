\name{cross.entropy}
\alias{cross.entropy}

\title{Cross-entropy criterion for snmf runs}
\description{
Return the cross-entropy criterion for runs of \code{snmf}cwith K ancestral populations.
The cross-entropy criterion is based on the prediction of masked genotypes to evaluate the fit of a model with K populations. The cross-entropy criterion helps choosing the number of ancestral populations or a best run for a fixed value of K. A smaller value of cross-entropy 
means a better run in terms of prediction capability.
The cross-entropy criterion is computed by the \code{\link{snmf}} function when the \code{entropy} Boolean option is TRUE.
}
\usage{
cross.entropy(object, K, run)
}

\arguments{
\item{object}{
A snmfProject object.
}
\item{K}{
The number of ancestral populations.
}
\item{run}{
A vector of run labels.
}
}
\value{
\item{res}{
A matrix containing the cross-entropy criterion for runs with K 
ancestral populations. 
}
}
\author{
    Eric Frichot
}
\seealso{
\code{\link{geno}}
\code{\link{snmf}}
\code{\link{G}}
\code{\link{Q}}
}
\examples{
### Example of analyses using snmf ###

# creation of a genotype file: genotypes.geno.
# The data contains 400 SNPs for 50 individuals.
data("tutorial")
write.geno(tutorial.R, "genotypes.geno")

################
# running snmf #
################

# Runs with K = 3 populations 
# cross-entropy is computed for 2 runs.
project = NULL
project = snmf("genotypes.geno", 
                K = 3, 
                entropy = TRUE, 
                repetitions = 2, 
                project = "new")

# get the cross-entropy for all runs for K = 3 
ce = cross.entropy(project, K = 3)

# get the cross-entropy for the 2nd run for K = 3
ce = cross.entropy(project, K = 3, run = 2)
}


\keyword{snmf}
