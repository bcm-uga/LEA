\name{ancestrymap}
\alias{ancestrymap}
\alias{example_ancestrymap}

\title{\code{ancestrymap} format description}
\description{
Description of the \code{ancestrymap} format. The \code{ancestrymap} format 
can be used as an input format for genotypic matrices in the functions 
\code{\link{pca}}, \code{\link{lfmm}} and \code{\link{snmf}}.
}
\details{
The \code{ancestrymap} format has one row for each genotype. Each row has 3 
columns: the 1st column is the SNP name, the 2nd column is the sample ID, 
the 3rd column is th number of alleles. Genotypes for a given SNP name 
are written in consecutive lines. The number of alleles can be the number 
of reference alleles or the number of derived alleles. Missing genotypes 
are encoded by the value 9.

Here is an example of a genotypic matrix using the 
ancestrymap format with 3 individuals and 4 SNPs:
\preformatted{
rs0000    SAMPLE0    1
rs0000    SAMPLE1    1
rs0000    SAMPLE2    2
rs1111    SAMPLE0    0
rs1111    SAMPLE1    1
rs1111    SAMPLE2    0
rs2222    SAMPLE0    0
rs2222    SAMPLE1    9
rs2222    SAMPLE2    1
rs3333    SAMPLE0    1
rs3333    SAMPLE1    2
rs3333    SAMPLE2    1
}
}

\author{
    Eric Frichot
}
\seealso{
\code{\link{ancestrymap2lfmm}}
\code{\link{ancestrymap2geno}}
\code{\link{geno}}
\code{\link{lfmm.data}}
\code{\link{ped}}
\code{\link{vcf}}
}


\keyword{format}
