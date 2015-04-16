\name{read.lfmm}
\alias{read.lfmm}

\title{Read files in the \code{\link{lfmm}} format}
\description{
Read a file in the \code{\link{lfmm}} format.
}
\usage{
read.lfmm(input.file)
}

\arguments{
\item{input.file}{
A character string containing a path to the input file,
a genotypic matrix in the \code{lfmm} format.
}
}
\value{
\item{R}{
A matrix containing the genotypes with one line per individual and one column 
per SNP.
}
}

\author{
    Eric Frichot
}
\seealso{
\code{\link{write.lfmm}}
\code{\link{lfmm.data}}
\code{\link{lfmm}}
\code{\link{geno2lfmm}}
\code{\link{lfmm2geno}}
\code{\link{ancestrymap2lfmm}}
\code{\link{ped2lfmm}}
}
\examples{
# tutorial contains a matrix of genotypes R with 1000 SNPs for 165 individuals.
# and a matrix with an environmental variable C.
data("tutorial")

# write R in a file called "genotypes.lfmm"
# Create file:    "genotypes.lfmm".
write.lfmm(tutorial.R,"genotypes.lfmm")

# read the file "genotypes.lfmm".
R = read.lfmm("genotypes.lfmm")
}


\keyword{read/write}

