\name{vcf}
\alias{vcf}
\alias{example_vcf}

\title{\code{vcf} format description}
\description{
Description of the \code{vcf} format. 
The \code{vcf} format can be used as an input format 
for genotypic matrices in the functions \code{\link{snmf}},
\code{\link{lfmm}}, and \code{\link{pca}}. 
}
\details{
The vcf format is described 
\href{http://www.1000genomes.org/wiki/Analysis/Variant\%20Call\%20Format/vcf-variant-call-format-version-41}{here}.

Here is an example of a genotypic matrix using the vcf format with 3 
individuals and 4 loci: 
\preformatted{
##fileformat=VCFv4.1 
##FORMAT=<ID=GM,Number=1,Type=Integer,Description="Genotype meta"> 
##INFO=<ID=VM,Number=1,Type=Integer,Description="Variant meta"> 
##INFO=<ID=SM,Number=1,Type=Integer,Description="SampleVariant meta"> 
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE0 SAMPLE1 SAMPLE2 
1 1001 rs0000 T C 999 . VM=1;SM=100 GT:GM 1/0:1 0/1:2 1/1:3 
1 1002 rs1111 G A 999 . VM=2;SM=101 GT:GM 0/0:6 0/1:7 0/0:8 
1 1003 notres G AA 999 . VM=3;SM=102 GT:GM 0/0:11 ./.:12 0/1:13 
1 1004 rs2222 G A 999 . VM=3;SM=102 GT:GM 0/0:11 . 1/0:13
1 1003 notres GA A 999 . VM=3;SM=102 GT:GM 0/0:11 ./.:12 0/1:13 
1 1005 rs3333 G A 999 . VM=3;SM=102 GT:GM 1/0:11 1/1:12 0/1:13 
}
}

\author{
    Eric Frichot
}
\seealso{
\code{\link{vcf2geno}}
\code{\link{vcf2lfmm}}
\code{\link{geno}}
\code{\link{lfmm}}
\code{\link{ped}}
\code{\link{ancestrymap}}
}


\keyword{format}
