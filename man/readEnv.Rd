\name{read.env}
\alias{read.env}

\title{Read environmental file in the \code{\link{env}}format}
\description{
Read a file in the \code{\link{env}} format.
}
\usage{
read.env(input.file)
}
\arguments{
\item{input.file}{
A character string containing a path to the input file,
an environmental data matrix in the \code{\link{env}} format.
}
}
\value{
\item{R}{
    A matrix containing the environmental variables with one line for each 
    individual and one column for each environmental variable.
}
}

\author{
    Eric Frichot
}
\seealso{
\code{\link{env}}
\code{\link{write.env}}
\code{\link{lfmm}}
}
\examples{
# Creation of an environmental matrix, C 
# containing 2 environmental variables for 3 individuals.
# C contains one line for each individual and one column for each variable.
C = matrix(runif(6), ncol=2, nrow=3)

# Write C in a file called "example.env".
# Create file:    "example.env".
write.env(C,"example.env")

# Read the file "example.env".
C = read.env("example.env")
}


\keyword{read/write}
