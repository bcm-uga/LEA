\name{Q}
\alias{Q}

\title{Admixture coefficients from a snmf run}
\description{
Return the \code{snmf} output matrix of admixture coefficients for the chosen 
run with K ancestral populations. For an example, see \code{\link{snmf}}.
}
\usage{
Q(object, K, run)
}

\arguments{
\item{object}{
A snmfProject object.
}
\item{K}{
The number of ancestral populations.
}
\item{run}{
A chosen run.
}
}
\value{
\item{res}{
A matrix containing the admixture coefficients for the chosen run with K 
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
\code{\link{cross.entropy}}
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
project = snmf("genotypes.geno", K = 3, repetitions = 2, project = "new")

# get the ancestry coefficients for the 2nd run for K = 3. 
res = Q(project, K = 3, run = 2)

# plot the 2nd run for K = 3 (ancestry coefficients).
barplot(t(Q(project, K = 3, run = 2)))
}


\keyword{snmf}
