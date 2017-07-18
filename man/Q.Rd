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
### Example of analysis using snmf ###

# Creation of the genotype file: genotypes.geno.
# The data contain 400 SNPs for 50 individuals.
data("tutorial")
write.geno(tutorial.R, "genotypes.geno")

################
# running snmf #
################

project.snmf <- snmf("genotypes.geno", 
                K = 3, 
                repetitions = 2, 
                project = "new")

# get the ancestry coefficients for the 2nd run for K = 3. 
Q.3 <- Q( project.snmf, K = 3, run = 2)

# cluster assignment for each individual
cluster <- apply( Q.3, 1, which.max)
table(cluster)
}


\keyword{snmf}
