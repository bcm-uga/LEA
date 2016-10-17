\name{ped}
\alias{ped}
\alias{example_ped}

\title{\code{ped} format description}
\description{
Description of the \code{ped} format. 
The \code{ped} format can be used as an input format 
for genotypic matrices in the functions \code{\link{snmf}},
\code{\link{lfmm}}, and \code{\link{pca}}. 
}
\details{
The \code{ped} format has one row for each individual. Each row contains 6
columns of information for each individual, plus two genotype columns for each
SNP. Each column must be separated by spaces or tabulations. The genotype 
format must be either 0ACGT or 01234, where 0 means missing genotype. The first
6 columns of the genotype file are: the 1st column is the family ID, the 2nd 
column is the sample ID, the 3rd and 4th columns are the sample IDs of parents,
the 5th column is the gender (male is 1, female is 2), the 6th column is the 
case/control status (1 is control, 2 is case), the quantitative trait value or
the population group label.

The ped format is described 
\href{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml\#ped}{here}. 

Here is an example with 3 individuals and 4 SNPs:
\preformatted{
1    SAMPLE0 0 0 2 2 1 2 3 3 1 1 2 1
2    SAMPLE1 0 0 1 2 2 1 1 3 0 4 1 1
3    SAMPLE2 0 0 2 1 2 2 3 3 1 4 1 2
}
}

\author{
    Eric Frichot
}
\seealso{
\code{\link{ped2lfmm}}
\code{\link{ped2geno}}
\code{\link{geno}}
\code{\link{lfmm.data}}
\code{\link{ancestrymap}}
\code{\link{vcf}}
}


\keyword{format}
