
vcf2lfmm <- function(    input.file,
            output.file = NULL,
            force = TRUE) 
{
    # test arguments and init
    # input file
        if(missing(input.file)) 
        stop("'input.file' argument is missing.")
    else if (!is.character(input.file))
        stop("'input.file' argument has to be of type character.")
    # check the extension
    test_extension(input.file, "vcf");
    # output file    
    if (!missing(output.file) && !is.character(output.file))
        stop("'output.file' argument has to be of type character.")
    else if (missing(output.file))
        output.file = setExtension(input.file, ".geno")
        # skip
        if (!force && file.exists(input.file) && file.exists(output.file)) {
#                print(cat("'", output.file, "' already exists!", sep=""))
                return(output.file)
    }

    # conversion
    tmp.file = vcf2geno(input.file);
    output.file = geno2lfmm(tmp.file);

    # create output 
    output.file;
}
