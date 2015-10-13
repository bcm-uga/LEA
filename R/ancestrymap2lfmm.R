ancestrymap2lfmm <- function(input.file = NULL, output.file = NULL, 
    force = TRUE)
{
    # test arguments and init
    # input file
        if(missing(input.file)) 
        stop("'input.file' argument is missing.")
    else if (!is.character(input.file))
        stop("'input.file' argument has to be of type character.")
    # check the extension
    test_extension(input.file, "ancestrymap");
    # output file    
    if (!missing(output.file) && !is.character(output.file))
        stop("'output.file' argument has to be of type character.")
        else if (missing(output.file)) 
        output.file = setExtension(input.file, ".lfmm")
        # skip
        if (!force && file.exists(input.file) && file.exists(output.file)) {
#               print(cat("'", output.file, "' already exists!", sep=""))
                return(output.file)
    }

    N = 0; M = 0;
    # run method
    .C("R_ancestrymap2lfmm", 
    as.character(input.file),
    as.character(output.file),
    as.integer(N),
    as.integer(M)
    );

    # create output 
    output.file;
}
