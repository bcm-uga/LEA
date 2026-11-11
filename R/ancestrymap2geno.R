
ancestrymap2geno <- function(input.file, output.file = NULL, force = TRUE)
{
    # test arguments and init
    # input file
    if(missing(input.file)) 
        stop("'input.file' argument is missing.")
    else if (!is.character(input.file))
        stop("'input.file' argument has to be of type character.")
    if (length(input.file) != 1L || is.na(input.file) || !nzchar(input.file))
        stop("'input.file' must be one non-missing, non-empty filename.")
    # check the extension
    test_extension(input.file, "ancestrymap");
    # output file    
    if (!missing(output.file) && !is.character(output.file))
        stop("'output.file' argument has to be of type character.")
    else if (missing(output.file)) 
        output.file = setExtension(input.file, ".geno")
    if (length(output.file) != 1L || is.na(output.file) || !nzchar(output.file))
        stop("'output.file' must be one non-missing, non-empty filename.")
    # skip
    if (!force && file.exists(input.file) && file.exists(output.file)) {
        return(output.file)
    }

    N = 0; M = 0;
    # run method
    .C("R_ancestrymap2geno", 
    as.character(input.file),
    as.character(output.file),
    as.integer(N),
    as.integer(M)
    );

    # create output 
    output.file;
}
