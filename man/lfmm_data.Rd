\name{lfmm.data}
\alias{lfmm.data}
\alias{example_lfmm}

\title{Input file for \code{\link{lfmm}}}
\description{
Description of the \code{lfmm} format.
The \code{lfmm} format can be used as an input format 
for genotypic matrices in the functions \code{\link{snmf}},
\code{\link{lfmm}}, \code{\link{lfmm2}}, and \code{\link{pca}}. 
}
\details{
The \code{lfmm} format has one row for each individual. Each row contains one 
value at each loci (separated by spaces or tabulations) corresponding to the 
number of alleles. The number of alleles corresponds to the number of reference
alleles or the number of derived alleles. Missing genotypes are encoded by the value -9 or the value 9. 

For the use of functions \code{\link{lfmm}} and \code{\link{lfmm2}} missing genotypes must be removed or imputed with the function \code{\link{impute}}. 

Here is an example of a genotypic matrix using the \code{lfmm} format with 3 individuals and 4 loci: 
\preformatted{
1 0 0 1
1 1 9 2
2 0 1 1
}
}
\author{
    Eric Frichot
}
\seealso{
\code{\link{lfmm}}
\code{\link{lfmm2}}
\code{\link{geno2lfmm}}
\code{\link{lfmm2geno}}
\code{\link{ancestrymap2lfmm}}
\code{\link{ped2lfmm}}
\code{\link{read.lfmm}}
\code{\link{write.lfmm}}
}


\keyword{format}
