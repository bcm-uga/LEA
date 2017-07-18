\name{snmf.pvalues}
\alias{snmf.pvalues}

\title{P-values for snmf population differentiation tests}
\description{
Returns a vector of p-values computed from an \code{snmf} run.
}
\usage{
snmf.pvalues (object, genomic.control, lambda, ploidy, entropy, K, run)
}

\arguments{
\item{object}{An snmfProject object.
}
\item{genomic.control}{
A Boolean value. If \code{TRUE}, the p-values are automatically calibrated using genomic control. If \code{FALSE}, the p-values are calculated
by rescaling the chi-squared test statistics using the \code{lambda} parameter.
}

\item{lambda}{A numeric value. The lambda value is used as an inflation factor to rescale the chi-squared statistics in the computation of p-values. This option requires that \code{genomic.control} = \code{FALSE}. The default value of \code{lambda} is equal to 1.0 (no rescaling).
}

\item{ploidy}{An integer value among 1 or 2. Tests are implemented for haploids and diploids (to be extended to polypoids).
}
\item{entropy}{A Boolean value. If \code{TRUE}, the run of minimum entropy is used for computing the p-values.
}
\item{K}{An integer value. The number of genetic clusters.
}
\item{run}{An integer for the run number used the computation of p-values (by default, the minimum entropy run).
}
}
\value{
\item{p.values}{A vector of p-values for each locus for the population differentiation test.}
\item{GIF}{The inflation factor value used in the test.}
}
\author{
    Olivier Francois
}
\seealso{
\code{\link{snmf}}
}

\references{
Martins, H., Caye, K., Luu, K., Blum, M. G. B., Francois, O. (2016). Identifying outlier loci in admixed and in continuous populations using ancestral population differentiation statistics. Molecular Ecology, 25(20), 5029-5042.
}


\examples{
### Example of analyses using snmf ###

data("tutorial")
# creation of a genotype file, "genotypes.lfmm".
# The data contain 400 SNPs for 50 individuals.
write.geno(tutorial.R, "genotypes.geno")


################
# snmf runs    #
################

# main options, K: the number of ancestral populations, 
#        entropy: cross-entropy criterion, 
#        CPU: the number of CPUs.

project.snmf = snmf("genotypes.geno", 
                    K = 4, 
                    entropy = TRUE,
                    ploidy = 2,
                    repetitions = 10,
                    project = "new")

# genome scan using adjusted p-values (genomic control method)

p = snmf.pvalues(project.snmf, entropy = TRUE, ploidy = 2, K = 4)
p$GIF

par(mfrow = c(2,1))
hist(p$pvalues, col = "orange")

plot(-log10(p$pvalues), pch = 19, col = "blue", cex = .7)

}


\keyword{snmf}
