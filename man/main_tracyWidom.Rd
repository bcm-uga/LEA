\name{tracy.widom}
\alias{tracy.widom}

\title{Tracy-Widom test for eigenvalues}
\description{
Perform tracy-widom tests on a set of eigenvalues to determine the number of 
significative eigenvalues and calculate the percentage of variance explained 
by each principal component. For an example, see \code{\link{pca}}. 
}
\usage{
tracy.widom (object)
}

\arguments{
\item{object}{
a pcaProject object.
}
}
\value{
\code{tracy.widom} returns a list containing the following components:
\item{eigenvalues}{
    The sorted input vector of eigenvalues (by descreasing order).
}
\item{twstats}{
    The vector of tracy-widom statistics.
}
\item{pvalues}{
    The vector of p-values associated with each eigenvalue. 
}
\item{effecn}{
    The vector of effective sizes.
}
\item{percentage}{
    The vector containing the percentage of variance 
    explained by each principal component.
}
}
\author{
    Eric Frichot
}
\references{
Tracy CA and Widom H. (1994). \emph{Level spacing distributions and the bessel 
kernel}. Commun Math Phys. 161 :289--309.
Patterson N, Price AL and Reich D. (2006). \emph{Population structure and 
eigenanalysis}. PLoS Genet. 2 :20.
}
\seealso{
\code{\link{pca}}
\code{\link{lfmm.data}}
\code{\link{lfmm}}
}
\examples{
# Creation of the genotype file "genotypes.lfmm"
# with 1000 SNPs for 165 individuals.
data("tutorial")
write.lfmm(tutorial.R,"genotypes.lfmm")

#################
# Perform a PCA #
#################

# run of PCA
# Available     options, K (the number of PCs calculated), 
#               center and scale. 
# Creation of   genotypes.pcaProject - the pcaProject object.
#               a directory genotypes.pca containing:
# Create files: genotypes.eigenvalues - eigenvalues,    
#               genotypes.eigenvectors - eigenvectors,
#               genotypes.sdev - standard deviations,
#               genotypes.projections - projections,
# Create a pcaProject object: pc.
pc = pca("genotypes.lfmm", scale = TRUE)

#############################
# Perform Tracy-Widom tests #
#############################

# Perfom Tracy-Widom tests on all eigenvalues.
# Create file:  genotypes.tracyWidom - tracy-widom test information, 
#               in the directory genotypes.pca/.
tw = tracy.widom(pc)

# Plot the percentage of variance explained by each component.
plot(tw$percentage)

# Display the p-values for the Tracy-Widom tests. 
tw$pvalues
}


\keyword{pca}
