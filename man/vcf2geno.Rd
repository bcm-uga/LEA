\name{vcf2geno}
\alias{vcf2geno}

\title{Convert from \code{\link{vcf}} to \code{\link{geno}} format}
\description{
A function that converts from the \code{\link{vcf}} format to the 
\code{\link{geno}} format. Note: This function may be obsolete. Conversion in accepted format such as \code{ped} can be obtained with the program \code{vcftools}.
}
\usage{
vcf2geno(input.file, output.file = NULL, force = TRUE)
}

\arguments{
\item{input.file}{
    A character string containing a path to the input file,
    a genotypic matrix in the \code{\link{vcf}} format. 
}
\item{output.file}{
    A character string containing a path to the output file,
    a genotypic matrix in the \code{\link{geno}} format. 
    By default, the name of the output file is the same name 
    as the input file with a .geno extension.
}
\item{force}{
    A boolean option. If FALSE, the input file is converted only if the output
    file does not exist. If TRUE, convert the file anyway.
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
\code{\link{vcf}}
\code{\link{geno}}
\code{\link{ancestrymap2lfmm}}
\code{\link{ancestrymap2geno}}
\code{\link{ped2lfmm}}
\code{\link{ped2geno}}
\code{\link{lfmm2geno}}
\code{\link{geno2lfmm}}
}
\examples{
# Creation of a file called "example.vcf"
# with 4 SNPs for 3 individuals.
data("example_vcf")
write.table(example_vcf,"example.vcf",col.names =
    c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
    "FORMAT", "SAMPLE0", "SAMPLE1", "SAMPLE2"),
    row.names = FALSE, quote = FALSE)

# Conversion    from the vcf format ("example.vcf") 
#               to the geno format ("example.geno").
# By default,   the name of the output file is the same name
#               as the input file with a .geno extension.
# Create files: "example.geno",
#               "example.vcfsnp" - SNP informations,
#               "example.removed" - removed lines. 
output = vcf2geno("example.vcf")

# Conversion    from the vcf format ("example.vcf")
#               to the geno format with the output file called "plop.geno".
# Create files: "plop.geno",
#               "plop.vcfsnp" - SNP informations,
#               "plop.removed" - removed lines. 
output = vcf2geno("example.vcf", "plop.geno")

# As force = false and the file "example.geno" already exists,
# nothing happens.
output = vcf2geno("example.vcf", force = FALSE)
}


\keyword{conversion}
