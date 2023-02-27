\name{genetic.gap}
\alias{genetic.gap}

\title{Genetic gap: genetic offset and genetic distance between environments.}

\description{The function returns estimates of the geometric genetic offset (genetic gap) computed from grids of new and predicted environments. The estimates are based on the covariance matrix of effect sizes obtained from an \code{lfmm2} model. The function takes as input the data that are used to adjust the LFMM, a matrix of environmental variables measured at new locations (new.env) or at the same locations as in the LFMM estimates (new.env = env is accepted), and a matrix of predicted environmental variables for the new locations (pred.env) in the same format as the new.env ones. 
}

\usage{
genetic.gap (input, env, new.env, pred.env, K, scale, candidate.loci)
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

\item{K}{An integer or a sequence of integers corresponding to the number of latent factors in the LFMM. The number of latent factors  could be estimated from the elbow point in the PCA screeplot for the genotype matrix. For a sequence of values, an average prediction will be returned.}

\item{scale}{A logical value indicating whether the environmental data are scaled or not. If \code{scale == TRUE}, then all environmental matrices are centered and scaled from the columwise mean and standard deviations of the \code{env} matrix. This option should be used only to evaluate the relative importance of environmental variables with the eigenvalues of the covariance matrix of effect sizes when the environmental data have different scales.}

\item{candidate.loci}{A vector specifying which loci (column label) in the genotype matrix are included in the computation of the genetic offset. The default value includes all loci.}
}

\value{
\item{offset}{A vector of genomic offset values computed for every sample location in \code{new.env} and \code{pred.env}. The genomic offset is the genetic gap defined in (Gain et al. 2023).
}
\item{distance}{A vector of environmental distance values computed for every sample location in \code{new.env} and \code{pred.env}. The distances to an estimate of the risk of nonadaptedness that includes correction for confounding factors and analyzes multiple predictors simultaneously (modified version of RONA). 
}
\item{eigenvalues}{Eigenvalues of the covariance matrix of LFMM effect sizes. They represent the relative importance of combinations of environmental variables described in \code{vectors} when the environmental data have similar scales. To be used with \code{scale == TRUE}.}

\item{vectors}{Eigenvectors of the covariance matrix of LFMM effect sizes representing combinations of environmental variablessorted by importance (\code{eigenvalues}).}
}

\author{
    Olivier Francois, Clement Gain
}

\seealso{
\code{\link{lfmm.data}}
\code{\link{lfmm2}}
}

\references{
Gain, C., et al. (2023). A quantitative theory for genomic offset statistics. bioRxiv, 10.1101/2023.01.02.522469.


Gain C, Francois O. (2021). LEA 3: Factor models in population genetics and ecological genomics with R. Molecular Ecology Resources. Molecular Ecology Resources 21 (8), 2738-2748. doi.org/10.1111/1755-0998.13366.
}

\examples{
### Example of genetic offset computation using lfmm2 ###

data("offset_example")

Y <- offset_example$geno
X <- offset_example$env
X.pred <- offset_example$env.pred

#PCA of the genotype data suggests k = 2 factors
plot(prcomp(Y), col = "blue")

## genetic gap

g.gap <- genetic.gap(input = Y, 
                           env = X, 
                           pred.env = X.pred,
                           K = 2)
                           
# return the values of the offset (genetic gap) for each sample location
round(g.gap$offset, digit = 3) 

# plot the squared root of the genetic gap vs Euclidean environmental distance
Delta = X - X.pred
dist.env = sqrt( rowSums(Delta^2) )  
plot(dist.env, sqrt(g.gap$offset),  cex = .6)


# plot RONA vs the genetic gap 
plot(g.gap$offset, g.gap$distance,  cex = .6)

# with scaled variables
g.gap.scaled <- genetic.gap(input = Y, 
                           env = X, 
                           pred.env = X.pred,
                           scale = TRUE,
                           K = 2)
                           
# Scaling does not change genetic gaps
plot(g.gap$offset, g.gap.scaled$offset,  cex = .6)                           
                           
# But scaling is useful for evaluating the relative importance of environmental variables
# Only two dimensions of the environmental space influence the genetic gap
barplot(g.gap.scaled$eigenvalues, col = "orange", xlab = "Axes", ylab = "Eigenvalues")

# The loadings for the first two variables indicate their relative contribution to local adaptation
g.gap.scaled$vectors[,1:2]

#rm(list = ls())
}

\keyword{lfmm2}


