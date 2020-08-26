\name{lfmm2.test}
\alias{lfmm2.test}

\title{P-values adjusted for latent factors computed by \code{lfmm2}.}

\description{The function returns a vector of p-values for association between loci and environmental variables adjusted for latent factors computed by \code{lfmm2}. It takes an object of class \code{lfmm2Class} with the data that were used to adjust the LFMM. 
}

\usage{
lfmm2.test (object, input, env, genomic.control, linear, family)
}

\arguments{
\item{object}{
An object of class \code{lfmm2Class}.
}

\item{input}{A genotypic matrix or a character string containing a path to the input file. The genotypic  matrix must be in the \code{\link{lfmm}{lfmm_format}} format without missing values (9 or \code{NA}). See \code{\link{impute}} for completion based on nonnegative matrix factorization and consider R packages for reading large matrices.
}

\item{env}{A matrix of environmental covariates or a character string containing a path to the environmental file. The environment matrix must be in the \code{\link{env}} format without missing values. Variables must be encoded as \code{numeric}.
}

\item{genomic.control}{
A logical value. If \code{TRUE}, the p-values are recalibrated by using genomic control after correction for confounding. 
}

\item{linear}{A logical value indicating whether linear or generalized linear models should be used to perform the association tests. If \code{FALSE}, \code{family} should be provided in the next argument.
}

\item{family}{a \code{family} for generalized linear models used in the association tests. The default is binomial(link = "logit")), which requires that y is between 0 and 1. 
}
}

\value{
\item{pvalues}{A matrix of p-values for each locus and each environmental variable.}
\item{zscores}{A matrix of z-scores for each locus and each environmental variable.}
\item{gif}{A vector of genomic inflation factors computed for each environmental variable.}
}

\author{
    Olivier Francois
}

\seealso{
\code{\link{lfmm.data}}
\code{\link{lfmm2}}
}

\references{
Caye K, Jumentier B, Lepeule J, Francois O. (2019). LFMM 2: fast and accurate inference of gene-environment associations in genome-wide studies. Molecular biology and evolution, 36(4), 852-860.
}

\examples{
### Example of analysis using lfmm2 ###

# Simulation with 10 target loci, with effect sizes ranging between -10 an 10 
# n = 100 individuals and L = 1000 loci

X <- as.matrix(rnorm(100)) # environmental variable
B <- rep(0, 1000) 
target <- sample(1:1000, 10) # target loci
B[target] <- runif(10, -10, +10) # effect sizes

# Creating hidden factors and loadings

U <- t(tcrossprod(as.matrix(c(-1,0.5,1.5)), X)) + matrix(rnorm(300), ncol = 3)
V <- matrix(rnorm(3000), ncol = 3)

# Simulating a binarized matrix containing haploid genotypes 
# Simulation performed with the generative LFMM

Y <- tcrossprod(as.matrix(X), B) + tcrossprod(U, V) + matrix(rnorm(100000, sd = .5), nrow = 100)
Y <- matrix(as.numeric(Y > 0), ncol = 1000)

######################################
# Fitting an LFMM with K = 3 factors #
######################################

mod2 <- lfmm2(input = Y, env = X, K = 3)

# Computing P-values and plotting their minus log10 values 
# Target loci are highlighted

pv <- lfmm2.test(object = mod2, input = Y, env = X, linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .4, pch = 19)
points(target, -log10(pv$pvalues[target]), col = "red")

}

\keyword{lfmm2}


