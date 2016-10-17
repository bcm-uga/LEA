read.lfmm <- function(input.file) {

        # test arguments
        if(missing(input.file))
                stop("'input.file' argument is missing.")
        else if (!is.character(input.file))
                stop("'input.file' argument has to be of type character.")
    # check extension 
    test_extension(input.file, "lfmm")

    return(as.matrix(read.table(input.file)))
}
