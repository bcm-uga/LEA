\name{read.geno}
\alias{read.geno}

\title{read a file in the \code{\link{geno}} format}
\description{
Read a file in the \code{\link{geno}} format.
}
\usage{
read.geno(input.file)
}

\arguments{
\item{input.file}{
A character string containing a path to the input file,
a genotypic matrix in the \code{geno} format.
}
}
\value{
\item{R}{
A matrix containing the genotypes with one line for each individual and one 
column for each SNP.
}
}
\author{
    Eric Frichot
}
\seealso{
\code{\link{write.geno}}
\code{\link{geno}}
\code{\link{snmf}}
\code{\link{geno2lfmm}}
\code{\link{lfmm2geno}}
\code{\link{ancestrymap2geno}}
\code{\link{ped2geno}}
\code{\link{vcf2geno}}
}
\examples{
# tutorial contains a matrix of genotypes R with 1000 SNPs for 165 individuals.
# and a matrix with an environmental variable C.
data("tutorial")

# Write R in a file called "genotypes.geno".
# Create file:    "genotypes.geno".
write.geno(tutorial.R,"genotypes.geno")

# Read the file "genotypes.geno".
R = read.geno("genotypes.geno")
}


\keyword{read/write}
