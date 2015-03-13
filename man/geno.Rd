\name{geno}
\alias{geno}
\alias{example_geno}

\title{Input file for \code{\link{snmf}}}
\description{
Description of the \code{geno} format.
The \code{geno} format can be used as an input format 
for genotypic matrices in the functions \code{\link{snmf}},
\code{\link{lfmm}}, and \code{\link{pca}}. 
}
\details{
The \code{geno} format has one row for each SNP.
Each row contains 1 character for each individual:
    0 means zero copy of the reference allele.
    1 means one copy of the reference allele.
    2 means two copies of the reference allele.
    9 means missing data.

Here is an example of a genotypic matrix using the \code{geno} format with 3 
individuals and 4 loci: 
\preformatted{
112
010
091
121
}
}

\author{
    Eric Frichot
}
\seealso{
\code{\link{geno2lfmm}}
\code{\link{lfmm2geno}}
\code{\link{ancestrymap2geno}}
\code{\link{ped2geno}}
\code{\link{vcf2geno}}
\code{\link{read.geno}}
\code{\link{write.geno}}
}


\keyword{format}
