\name{cross.entropy.estimation}
\alias{cross.entropy.estimation}

\title{compute the cross-entropy criterion}
\description{
Calculate the cross-entropy criterion. This is an internal function, 
automatically called by \code{\link{snmf}}.
The cross-entropy criterion is a value based on the prediction of masked 
genotypes to evaluate the error of ancestry estimation. The criterion will help
to choose the best number of ancestral population (K) and the best run among a
set of runs in \code{\link{snmf}}. A smaller value of cross-entropy means a 
better run in terms of prediction capacity. 
The cross.entropy.estimation function displays the cross-entropy criterion 
estimated on all data and on masked data based on the input file, the masked 
data file (created by \code{\link{create.dataset}}, the estimation of the 
ancestry coefficients Q and the estimation of ancestral genotypic frequencies, 
G (calculated by \code{\link{snmf}}).  
The cross-entropy estimation for all data is always lower than the 
cross-entropy estimation for masked data. The cross-entropy estimation useful
to compare runs is the cross-entropy estimation for masked data.
The cross-entropy criterion can also be automatically calculated by the 
\code{\link{snmf}} function with the \code{entropy} option. 
}
\usage{
cross.entropy.estimation (input.file, K, masked.file, Q.file, G.file, 
    ploidy = 2)
}

\arguments{
\item{input.file}{
    A character string containing a path to the input file without masked 
    genotypes, a genotypic matrix in the \code{\link{geno}} format.
}
\item{K}{
    An integer corresponding to the number of ancestral populations.
}
\item{masked.file}{
    A character string containing a path to the input file with masked 
    genotypes, a genotypic matrix in the \code{\link{geno}} format. This file 
    can be generated with the function, \code{\link{create.dataset}}).
    By default, the name of the masked data file is the same name as the input
    file with a \code{_I.geno} extension.
}
\item{Q.file}{
    A character string containing a path to the input ancestry coefficient 
    matrix Q. 
    By default, the name of this file is the same name as the input file with 
    a \code{K.Q} extension.
}
\item{G.file}{
    A character string containing a path to the input ancestral genotype 
    frequency matrix G. By default, the name of this file is the same name as 
    the input file with a \code{K.G} extension (\code{input_file.K.G}).
}
\item{ploidy}{
    1 if haploid, 2 if diploid, n if n-ploid.
}
}
\value{
\code{cross.entropy.estimation} returns a list containing the following 
components:
\item{masked.ce}{
    The value of the cross-entropy criterion of the masked genotypes.
}
\item{all.ce}{
    The value of the cross-entropy criterion of all the genotypes.
}
}
\author{
    Eric Frichot
}
\seealso{
\code{\link{geno}}
\code{\link{create.dataset}}
\code{\link{snmf}}
}
\references{
Frichot E, Mathieu F, Trouillon T, Bouchard G, Francois O. (2014). \emph{Fast 
and Efficient Estimation of Individual Ancestry Coefficients}. Genetics, 194(4)
: 973--983.
}
\examples{

# Creation of tuto.geno
# A file containing 400 SNPs for 50 individuals.
data("tutorial")
write.geno(tutorial.R,"genotypes.geno")

# The following command are equivalent with 
# project = snmf("genotypes.geno", entropy = TRUE, K = 3)
# cross.entropy(project)

# Creation      of the masked data file
# Create file:  "genotypes_I.geno"
output = create.dataset("genotypes.geno")

# run of snmf with genotypes_I.geno and K = 3
project = snmf("genotypes_I.geno", K = 3, project = "new")

# calculate the cross-entropy
res = cross.entropy.estimation("genotypes.geno", K = 3, "genotypes_I.geno",
    "./genotypes.snmf/K3/run1/genotypes_r1.3.Q", 
    "./genotypes.snmf/K3/run1/genotypes_r1.3.G")

# get the result
res$masked.ce
res$all.ce
}
