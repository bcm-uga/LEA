\name{G}
\alias{G}

\title{Ancestral allele frequencies from a snmf run}
\description{
Return the \code{snmf} output matrix of ancestral allele frequency matrix for 
the chosen run with K ancestral populations. For an example, see 
\code{\link{snmf}}.
}
\usage{
G(object, K, run)
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
A matrix containing the ancestral allele frequencies for a run with 
K ancestral populations. 
}
}
\author{
    Eric Frichot
}
\seealso{
\code{\link{geno}}
\code{\link{snmf}}
\code{\link{Q}}
\code{\link{cross.entropy}}
}
\examples{
### Example of analyses using snmf ###

# creation of a genotype file: genotypes.geno.
# The data contain 400 SNPs for 50 individuals.
data("tutorial")
write.geno(tutorial.R, "genotypes.geno")

################
# running snmf #
################


# Two runs for K = 1 to 5 
project.snmf = snmf("genotypes.geno", 
                K = 3, 
                repetitions = 2, 
                project = "new")

# get the ancestral genotype frequency matrix, G, for the 2nd run for K = 3. 
freq = G(project.snmf, K = 3, run = 2)
}


\keyword{snmf}
