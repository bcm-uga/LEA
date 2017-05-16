\name{barplot.Qmatrix}
\alias{barplot.Qmatrix}
\title{Barplot representation of a Q-matrix}
\usage{
\method{barplot}{Qmatrix}(height, sort.by.Q = TRUE, lab = FALSE, ...)
}
\arguments{
\item{height}{an object of class \code{Qmatrix} containing a matrix of ancestry coefficients computed from \code{snmf} or converted from another format.}

\item{sort.by.Q}{a Boolean value indicating whether individuals should be sorted by ancestry or not.}

\item{lab}{a list of individual labels.}

\item{...}{other parameters of the function \code{\link{barplot.default}}.}
}
\value{
A permutation of individual labels used in the sort.by.Q option (order). Displays the
Q matrix.
}
\description{
This function displays a barplot representation of the
ancestry coefficient matrix. It includes a sort-by-Q option.
}
\examples{
library(LEA)
# creation of a genotype file: genotypes.geno.
# 400 SNPs for 50 individuals.

data("tutorial")
write.geno(tutorial.R, "genotypes.geno")

################
# running snmf #
################

project = snmf("genotypes.geno", K = 4, entropy = TRUE, repetitions = 10,
    project = "new")

# get the cross-entropy value for each run 
ce = cross.entropy(project, K = 4)

# select the run with the lowest cross-entropy value
best = which.min(ce)

# plot the best run for K = 4 (ancestry coefficients).
Q.matrix <- as.qmatrix(Q(project, K = 4, run = best))
my.colors <- c("tomato", "lightblue", "olivedrab", "gold")

barplot(Q.matrix, border = NA, space = 0, col = my.colors, xlab = "Individuals",
        ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
        
axis(1, at = 1:nrow(Q.matrix), labels = bp$order, las = 3, cex.axis = .4)
}
\author{
Olivier François
}
\seealso{
\code{ \code{\link{as.qmatrix}} \code{\link{snmf}} }
}

