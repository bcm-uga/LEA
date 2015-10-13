\name{create.dataset}
\alias{create.dataset}

\title{create a data set with masked data}
\description{

\code{\link{create.dataset}} creates a data set with a given percentage of 
masked data from the original data set. It is used to calculate the 
\code{\link{cross.entropy}} criterion.
}
\usage{
create.dataset (input.file, output.file, seed = -1, percentage = 0.05)
}
\details{
This is an internal function, automatically called by \code{\link{snmf}} with 
the \code{entropy} option.
}

\arguments{
\item{input.file}{
    A character string containing a path to the input file,
    a genotypic matrix in the \code{\link{geno}} format.
}
\item{output.file}{
    A character string containing a path to the output file,
    a genotypic matrix in the \code{\link{geno}} format. 
    The output file is the input file with masked genotypes. 
    By default, the name of the output file is the same name as 
    the input file with a _I.geno extension.
}
\item{seed}{
    A seed to initialize the random number generator. 
    By default, the seed is randomly chosen.
}
\item{percentage}{
    A numeric value between 0 and 1 containing the percentage of masked 
    genotypes.  
}
}
\value{
\item{output.file}{
    A character string containing a path to the output file,
    a genotypic matrix in the \code{\link{geno}} format. 
}
}
\author{
    Eric Frichot
}
\seealso{
\code{\link{geno}}
\code{\link{snmf}}
\code{\link{cross.entropy}}
}
\examples{
# Creation of tuto.geno
# A file containing 400 SNPs for 50 individuals.
data("tutorial")
write.geno(tutorial.R,"genotypes.geno")

# Creation      of the masked data file
# Create file:  "genotypes_I.geno"
output = create.dataset("genotypes.geno")
}

