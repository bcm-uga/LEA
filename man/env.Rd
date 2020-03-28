\name{env}
\alias{env}

\title{Environmental input file format for \code{\link{lfmm}}}
\description{

Description of the \code{env} format. The \code{env} format 
can be used as an input format for the environmental variables 
in the \code{\link{lfmm}} function. 
}
\details{
The \code{env} format has one row for each individual. Each row contains one value for each environmental variable (separated by spaces or tabulations).

Here is an example of an environmental file using the \code{env} format with 3 individuals and 2 variables:
\preformatted{
0.252477 0.95250639
0.216618 0.10902647
-0.47509 0.07626694
}
}
\author{
    Eric Frichot
}
\seealso{
\code{\link{lfmm}}
\code{\link{lfmm2}}
\code{\link{read.env}}
\code{\link{write.env}}
}


\keyword{format}
