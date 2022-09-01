\name{genetic.gap}
\alias{genetic.gap}

\title{Genetic gap: genetic offset and genetic distance between environments.}

\description{The function returns genetic gap estimates computed from grids of new and predicted environments. The estimates are based on the covariance matrix of effect sizes obtained from an \code{lfmm2} model. It takes as input the data that are used to adjust the LFMM, a matrix of environmental variables measured at new locations (new.env) or at the same locations as in the LFMM estimates (new.env = env is accepted), and a matrix of predicted environmental variables for the new locations (pred.env) in the same format as the new.env ones. 
}

\usage{
genetic.gap (input, env, new.env, pred.env, K, candidate.loci)
}

\arguments{


\item{input}{A genotypic matrix or a character string containing a path to the input file. The genotypic  matrix must be in the \code{\link{lfmm}} format without missing values (9 or \code{NA}). See \code{\link{impute}} for completion based on nonnegative matrix factorization. Also consider R packages for reading large matrices.
}

\item{env}{A matrix of environmental covariates or a character string containing a path to the environmental file. The environmental matrix must be in the \code{\link{env}} format without missing values. The variables must be encoded as \code{numeric} and sampled at the same locations as for the genotype matrix.
}

\item{new.env}{A matrix of new environmental covariates or a character string containing a path to the new environmental data. The data are environmental covariates sampled at locations that can differ from those used in the estimation of the LFMM (env). By default, the matrix provided as the \code{env} argument is used.  The new environmental matrix must be in the \code{\link{env}} format without missing values. The variables must be encoded as \code{numeric} 
}

\item{pred.env}{A matrix of predicted (new) environmental covariates or a character string containing a path to the predicted environmental data file. The predicted environmental matrix must be in the \code{\link{env}} format without missing values, and of same dimension as the \code{new.env} matrix. All variables must be encoded as \code{numeric} and sampled at the same locations as for the \code{new.env} matrix. Predicted environmental covariates typically result from bioclimatic models (eg, worldclim). 
}

\item{K}{An integer corresponding to the number of latent factors in the LFMM. It could be estimated from the elbow point in the PCA screeplot for the \code{geno} (genotype) matrix.}

\item{candidate.loci}{A vector specifying which loci (column label) in the genotype matrix are included in the computation of the genetic offset (default value: all positions).}
}

\value{
\item{offset}{A vector of genetic gap values computed for every sample location in \code{new.env} and \code{pred.env}.
}
\item{distance}{A vector of environmental distance values computed for every sample location in \code{new.env} and \code{pred.env}. The distances correspond to the squared root of the offset measures. 
}
\item{eigenvalues}{Eigenvalues of the covariance matrix of LFMM effect sizes.}

\item{vectors}{Eigenvectors of the covariance matrix of LFMM effect sizes.}
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
### Example of genetic gap computation using lfmm2 ###

data("offset_example")

Y <- offset_example$geno
X <- offset_example$env
X.pred <- offset_example$env.pred

#PCA of the genotype data suggests k = 2 factors
plot(prcomp(Y), col = "blue")

g.gap <- genetic.gap(input = Y, 
                           env = X, 
                           pred.env = X.pred,
                           K = 2)
                           
# return the squared root of the genetic gap for each sample location
round(g.gap$distance, digit = 3) 

# plot the squared root of the genetic gap vs Euclidean environmental distance
dist.env = sqrt(diag((X - X.pred) %*% t(X - X.pred)))
plot(dist.env, g.gap$distance,  cex = .6)

#rm(list = ls())
}

\keyword{lfmm2}


