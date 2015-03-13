\name{geno2lfmm}
\alias{geno2lfmm}

\title{Convert from \code{\link{geno}} to \code{\link{lfmm}} format}
\description{
A function that converts from the \code{\link{geno}} format to the 
\code{\link{lfmm}} format.
}
\usage{
geno2lfmm(input.file, output.file = NULL, force = TRUE)
}

\arguments{
\item{input.file}{
    A character string containing a path to the input file,
    a genotypic matrix in the \code{\link{geno}} format. 
}
\item{output.file}{
    A character string containing a path to the output file,
    a genotypic matrix in the \code{\link{lfmm}} format. 
    By default, the name of the output file is the same name 
    as the input file with a .lfmm extension.
}
\item{force}{
    A boolean option. If FALSE, the input file is converted only if the output
    file does not exist. If TRUE, convert the file anyway.
}
}
\value{
\item{output.file}{
    A character string containing a path to the output file,
    a genotypic matrix in the \code{\link{lfmm}} format. 
}
}
\author{
    Eric Frichot
}
\seealso{
\code{\link{lfmm.data}}
\code{\link{geno}}
\code{\link{ancestrymap2lfmm}}
\code{\link{ancestrymap2geno}}
\code{\link{ped2lfmm}}
\code{\link{ped2geno}}
\code{\link{vcf2geno}}
\code{\link{lfmm2geno}}
\code{\link{read.geno}}
\code{\link{write.geno}}
}
\examples{
# Creation of a file called "genotypes.geno" in the working directory
# with 400 SNPs for 50 individuals.
data("tutorial")
write.geno(tutorial.R, "genotypes.geno")

# Conversion    from the geno format ("genotypes.geno") 
#               to the lfmm format ("genotypes.lfmm").
# By default,   the name of the output file is the same name 
#               as the input file with a .lfmm extension.
# Create file:  "genotypes.lfmm".
output = geno2lfmm("genotypes.geno")

# Conversion    from the geno format ("genotypes.geno") 
#               to the lfmm format with the output file called "plop.lfmm".
# Create file:  "plop.lfmm".
output = geno2lfmm("genotypes.geno", "plop.lfmm")

# As force = false and the file "genotypes.lfmm" already exists,
# nothing happens.
output = geno2lfmm("genotypes.geno", force = FALSE)
}


\keyword{conversion}
