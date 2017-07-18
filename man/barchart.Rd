\name{barchart}
\alias{barchart}
\title{Bar plot representation of an snmf Q-matrix}
\description{
This function displays a bar plot/bar chart representation of the Q-matrix computed from an \code{snmf} run. The function can use a sort by Q option. See \code{\link{snmf}}.
}

\usage{
barchart (object, K, run, sort.by.Q = TRUE, lab = FALSE, ...)
}
\arguments{
\item{object}{an snmfProject object.
}

\item{K}{an integer value corresponding to number of ancestral populations.
}

\item{run}{an integer value. Usually the run number that minimizes the cross-entropy criterion.
}

\item{sort.by.Q}{a Boolean value indicating whether individuals should be sorted by their 
ancestry or not.
}

\item{lab}{a list of individual labels.
}

\item{...}{other parameters of the function \code{\link{barplot.default}}.}
}
\value{
A permutation of individual labels used in the sort.by.Q option (order). Displays the
Q matrix.
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

# plot the ancestry coefficients for the best run and K = 4 

my.colors <- c("tomato", "lightblue", "olivedrab", "gold")
barchart(project, K = 4, run = best, border = NA, space = 0, col = my.colors, xlab = "Individuals",
        ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
        
axis(1, at = 1:nrow(Q.matrix), labels = bp$order, las = 3, cex.axis = .4)
}
\author{
Olivier François
}
\seealso{
\code{ \code{\link{snmf}} }
}

