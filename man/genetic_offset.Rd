\name{genetic.offset}
\alias{genetic.offset}

\title{Genetic offset for populations under new environments \code{lfmm2}.}
\description{The function returns a vector of genetic offset estimates computed from specified population labels and new environments using predictions of the \code{lfmm2} fitted model. It takes an object of class \code{lfmm2Class} with the data that were used to adjust the LFMM, plus new environmental variables in the same format as the original ones. 
}
\usage{
genetic.offset (object, input, env, new.env, pop.labels)
}

\arguments{
\item{object}{
An object of class \code{lfmm2Class}.
}

\item{input}{A genotypic matrix or a character string containing a path to the input file. The genotypic  matrix must be in the \code{\link{lfmm}{lfmm_format}} format without missing values (9 or \code{NA}). See \code{\link{impute}} for completion based on nonnegative matrix factorization and consider R packages for reading large matrices.
}

\item{env}{A matrix of environmental covariates or a character string containing a path to the environmental file. The environment matrix must be in the \code{\link{env}} format without missing values. All variables must be encoded as \code{numeric}.
}

\item{env}{A matrix of new environmental covariates or a character string containing a path to the new environmental data file. The new environment matrix must be in the \code{\link{env}} format without missing values, and of same dimension as the \code{env} matrix. All variables must be encoded as \code{numeric}.
}

\item{pop.labels}{A numeric or character vector providing a population label for each individual in the response matrix. 
}
}
\value{
\item{offset}{A vector of genetic offset values computed for each population.}
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
### Example of offset prediction using lfmm2 ###

# Simulation with 10 target loci, with effect sizes ranging between -10 an 10 
# n = 100 individuals and L = 1000 loci

X <- as.matrix(rnorm(100)) # environmental variable
B <- rep(0, 1000) 
target <- sample(1:1000, 100) # target loci
B[target] <- runif(10, -10, +10) # effect sizes

# Creating hidden factors and loadings

U <- t(tcrossprod(as.matrix(c(-1,0.5,1.5)), X)) + matrix(rnorm(300), ncol = 3)
V <- matrix(rnorm(3000), ncol = 3)

# Simulating a binarized matrix containing haploid genotypes 
# Simulation performed with a generative LFMM

Y <- tcrossprod(as.matrix(X), B) + tcrossprod(U, V) + matrix(rnorm(100000, sd = .5), nrow = 100)
Y <- matrix(as.numeric(Y > 0), ncol = 1000)

######################################
# Fitting an LFMM with K = 3 factors #
######################################

mod2 <- lfmm2(input = Y, env = X, K = 3)

# Computing genetic offsets for 2 populations, defined from latent factor 1
pop <-  1 + (U[,1] > 0)

g.offset <- genetic.offset(object = mod2, input = Y, 
                           env = X, new.env = X + 2, 
                           pop.labels = pop)
g.offset 

#rm(list = ls())
}
\keyword{lfmm2}


