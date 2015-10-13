\name{ancestrymap2lfmm}
\alias{ancestrymap2lfmm}

\title{Convert from \code{\link{ancestrymap}} to \code{\link{lfmm}} format}
\description{
A function that converts from the \code{\link{ancestrymap}} format to the 
\code{\link{lfmm}} format.
}
\usage{
ancestrymap2lfmm(input.file, output.file = NULL, force = TRUE)
}

\arguments{
\item{input.file}{
A character string containing a path to the input file, 
a genotypic matrix in the \code{\link{ancestrymap}} format. 
}
\item{output.file}{
A character string containing a path to the output file,
a genotypic matric in the \code{\link{lfmm}} format. 
By default, the name of the output file is the same name 
as the input file with a .lfmm extension.
}
\item{force}{
    A boolean option. If FALSE, the input file is converted only if the output
    file does not exist. If TRUE, convert the file anyway.
}
}
\value{
\item{output.file}{
A character string containing a path to the output file,
a genotypic matric in the \code{\link{lfmm}} format. 
}
}
\author{
    Eric Frichot
}

\seealso{
\code{\link{ancestrymap}}
\code{\link{lfmm.data}}
\code{\link{ancestrymap2geno}}
\code{\link{geno2lfmm}}
\code{\link{ped2lfmm}}
\code{\link{ped2geno}}
\code{\link{vcf2geno}}
\code{\link{lfmm2geno}}
}


\examples{
# Creation of a file called "example.ancestrymap"
# containing 4 SNPs for 3 individuals.
data("example_ancestrymap")
write.table(example_ancestrymap,"example.ancestrymap",
col.names = FALSE, row.names = FALSE, quote = FALSE)

# Conversion    from the ancestrymap format ("example.ancestrymap") 
#               to the lfmm format ("example.lfmm").
# By default,   the name of the output file is the same name
#               as the input file with a .lfmm extension.
# Create file:  "example.lfmm".
output = ancestrymap2lfmm("example.ancestrymap")

# Conversion    from the ancestrymap format (example.ancestrymap) 
#               to the geno format with the output file called plop.lfmm.
# Create file:  "plop.lfmm".
output = ancestrymap2lfmm("example.ancestrymap", "plop.lfmm")

# As force = false and the file "example.lfmm" already exists,
# nothing happens.
output = ancestrymap2lfmm("example.ancestrymap", force = FALSE)
}



\keyword{conversion}
