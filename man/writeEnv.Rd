\name{write.env}
\alias{write.env}

\title{Write files in the \code{\link{env}} format}
\description{
Write a file in the \code{\link{env}} format.
}
\usage{
    write.env(R, output.file)
}

\arguments{
\item{R}{
    A matrix containing the environmental variables with one line for each 
    individual and one column for each environmental variable. The missing 
    genotypes have to be encoded with the value 9.
}
\item{output.file}{
    A character string containing a path to the output file,
    an environmental data matrix in the \code{env} formt.
}
}
\value{
\item{output.file}{
    A character string containing a path to the output file,
    an environmental data matrix in the \code{env} formt.
}
}

\author{
    Eric Frichot
}
\seealso{
\code{\link{read.env}}
\code{\link{env}}
\code{\link{lfmm}}
}
\examples{
# Creation of an environmental matrix C 
# containing 2 environmental variables for 3 individuals.
# C contains one line for each individual and one column for each variable.
C = matrix(runif(6), ncol=2, nrow=3)

# Write C in a file called "tuto.env".
# Create file:    "tuto.env".
write.env(C,"tuto.env")

# Read the file "tuto.env".
C = read.env("tuto.env")
}


\keyword{read/write}
