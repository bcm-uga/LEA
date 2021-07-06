\name{impute}
\alias{impute}

\title{Impute missing genotypes using an snmf object}
\description{
Impute missing genotypes in a genotype file (.lfmm) by using ancestry and genotype frequency estimates from an \code{snmf} run. The function generates a new \code{lfmm} file. See \code{\link{lfmm}} and \code{\link{lfmm2}}.
}
\usage{
impute (object, input.file, method, K, run) 
}

\arguments{
\item{object}{
An snmfProject object.
}
\item{input.file}{A path (character string) to an input file in lfmm format with missing genotypes. The same input data must be used when generating the snmf object. 
}
\item{method}{A character string: "random" or "mode". With "random", imputation is performed by using the genotype probabilities. With "mode", the most likely genotype is used for matrix completion. 
}
\item{K}{An integer value. The number of ancestral populations.
}

\item{run}{An integer value. A particular run used for imputation (usually the run number that minimizes the cross entropy criterion).
}
}
\value{
\item{NULL}{The function writes the imputed genotypes in an output file having the "_imputed.lfmm" suffix.
}
}
\author{
  Olivier Francois
}
\seealso{
\code{\link{snmf}}
\code{\link{lfmm}}
\code{\link{lfmm2}}
}

\references{
Gain C, Francois O (2021). LEA 3: Factor models in population genetics and ecological genomics with R. Molecular Ecology Resources, 
doi.org/10.1111/1755-0998.13366. 
}


\examples{
### Example of analysis ###

data("tutorial")
# creation of a genotype file with missing genotypes
# The data contain 400 SNPs for 50 individuals.

dat = as.numeric(tutorial.R)
dat[sample(1:length(dat), 100)] <-  9
dat <- matrix(dat, nrow = 50, ncol = 400 )
write.lfmm(dat, "genotypes.lfmm")

################
# running snmf #
################

project.snmf = snmf("genotypes.lfmm", K = 4, 
        entropy = TRUE, repetitions = 10,
        project = "new")
        
# select the run with the lowest cross-entropy value
best = which.min(cross.entropy(project.snmf, K = 4))

# Impute the missing genotypes
impute(project.snmf, "genotypes.lfmm", method = 'mode', K = 4, run = best)

# Compare with truth
# Proportion of correct imputation results:
mean( tutorial.R[dat == 9] == read.lfmm("genotypes.lfmm_imputed.lfmm")[dat == 9] )

}


\keyword{snmf}
