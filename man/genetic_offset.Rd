\name{genetic.offset}
\alias{genetic.offset}

\title{Population genetic offset under new environments.}

\description{The function returns genetic offset estimates computed from user-specified population labels and new environments based on predictions of an \code{lfmm2} model. It takes as input the data that are used to adjust the LFMM, and a matrix of new (predicted) environmental variables in the same format as the original ones. 
}

\usage{
genetic.offset (input, env, pred.env, pop.labels, K, pca, candidate.loci)
}

\arguments{


\item{input}{A genotypic matrix or a character string containing a path to the input file. The genotypic  matrix must be in the \code{\link{lfmm}} format without missing values (9 or \code{NA}). See \code{\link{impute}} for completion based on nonnegative matrix factorization and consider R packages for reading large matrices.
}

\item{env}{A matrix of environmental covariates or a character string containing a path to the environmental file. The environment matrix must be in the \code{\link{env}} format without missing values. The variables must be encoded as \code{numeric} and sampled at the same locations as for the genotype matrix.
}

\item{pred.env}{A matrix of predicted (new) environmental covariates or a character string containing a path to the predicted environmental data file. The predicted environmental matrix must be in the \code{\link{env}} format without missing values, and of same dimension as the \code{env} matrix. All variables must be encoded as \code{numeric}. Examples of predicted environmental covariates could result from bioclimatic models (eg, worldclim).
}

\item{pop.labels}{A numeric or character vector providing population labels for all rows (individuals) of the response matrix.}

\item{K}{An integer corresponding to the number of latent factors in the LFMM. The \code{NULL} default takes K equal to the number of samples minus the number of environmental covariates.}

\item{pca}{A boolean value indicating whether genetic offsets are computed from a PCA approximation or not (default value: \code{FALSE}). 
}

\item{candidate.loci}{A vector specifying which loci (column label) in the genotype matrix are included in the computation of the genetic offset (default value: all positions).}

}

\value{
\item{offset}{A matrix or vector (depending on the \code{pca} parameter) of genetic offset values computed for every population in \code{pop.labels}.}
}

\author{
    Olivier Francois, Clement Gain
}

\seealso{
\code{\link{lfmm.data}}
\code{\link{lfmm2}}
}

\references{
Gain C, Francois O. (2021). LEA 3: Factor models in population genetics and ecological genomics with R. Molecular Ecology Resources. Molecular Ecology Resources 21 (8), 2738-2748. doi.org/10.1111/1755-0998.13366.
}

\examples{
### Example of offset prediction using lfmm2 ###

# Simulation with 100 target loci
# Effect sizes ranging between -10 an 10 
# n = 100 individuals and L = 1000 loci

X <- as.matrix(rnorm(100)) # environmental variable
B <- rep(0, 1000) 
target <- sample(1:1000, 100) # target loci
B[target] <- runif(100, -10, +10) # effect sizes

# Creating hidden factors and their loadings

U <- t(tcrossprod(as.matrix(c(-1.25,0.5,1.25)), X)) + 
      matrix(rnorm(300), ncol = 3)
      
V <-  matrix(rnorm(3000), ncol = 3)

# Simulating a binarized matrix of haploid genotypes 
# Simulation performed with a generative LFMM

Y <- tcrossprod(as.matrix(X), B) + tcrossprod(U, V) + 
      matrix(rnorm(100000, sd = .5), nrow = 100)
      
Y <-  matrix(as.numeric(Y > 0), ncol = 1000)

######################################
# Fitting an LFMM with K = 98 factors#
######################################


# Computing genetic offset statistics for 2 populations 
# defined from PCA

pop <- 1 + (prcomp(Y, scale = TRUE)$x[,1] > 0)

g.offset <- genetic.offset(input = Y, 
                           env = X, pred.env = 2*X + 10, 
                           pop.labels = pop)
round(g.offset, digit = 3) 



# Computing genetic offset statistics for 2 populations, 
# defined from ANOVA and  truly targeted loci

g.offset <- genetic.offset(input = Y, 
                           env = X, pred.env = 2*X + 10, 
                           pop.labels = pop,
                           pca = FALSE, 
                           candidate.loci = target)
                           
round(g.offset, digit = 3) 

#rm(list = ls())
}

\keyword{lfmm2}


