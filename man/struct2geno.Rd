\name{struct2geno}
\alias{struct2geno}

\title{Conversion from the STRUCTURE format to the geno format.}
\description{
The function converts a multiallelic genotype file in the STRUCTURE format into a file in the 'geno' for \code{\link{snmf}} and the 'lfmm' format for \code{\link{lfmm}}.
}
\usage{
struct2geno (input.file, ploidy, FORMAT, extra.row, extra.column)
}
\arguments{
\item{input.file}{A character string. A path to a STRUCTURE or a TESS input file of multiallelic markers (eg, microsatellites) for haploid or diploid individuals. Missing data must be encoded as "-9" or as any negative value. Individual genotypes are encoded using either one or two rows of data. 
}
\item{genomic.control}{
A Boolean value. If \code{TRUE}, the p-values are automatically calibrated using genomic control. If \code{FALSE}, the p-values are calculated
by rescaling the chi-squared test statistics using the \code{lambda} parameter.
}
\item{ploidy}{An integer value (1 or 2). Value 2 for diploids and 1 for haploids.
}
\item{FORMAT}{An integer value equal to 1 for markers encoded using one row of data for each individual, and 2 for markers encoded using two rows of data for each individual.
}
\item{extra.row}{An integer value indicating the number of extra rows in the header of the input file (eg, marker ids).
}
\item{extra.column}{an integer value indicating the number of extra columns in the input file. Extra columns can include individual ids, pop ids, geographic coordinates, etc.
}
}
\value{
NULL. Output files in the 'geno' and the 'lfmm' format record individual genotypes for each allele at each marker.
}
\author{
    Olivier François
}
\seealso{
\code{\link{lfmm.data}}
\code{\link{geno}}
\code{\link{lfmm}}
\code{\link{snmf}}
}
\examples{
### Example of conversion from a STRUCTURE format ###
### Artificial data with 10 diploid individuals and 10 STR markers
### Input file: 'dat.str'

library(LEA)
dat.str  <- matrix(sample(c(100:120,-9), 200, replace = TRUE), 
                   nrow = 10, ncol = 20)
write.table(dat.str, 
            file = "dat.str", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)

### Conversion 
struct2geno("dat.str", ploidy = 2, FORMAT = 1)

### snmf run and barplot
s  <- snmf("dat.str.geno", K = 2, project = "new")
barplot(as.qmatrix(Q(s)), xlab = "Individuals")
}