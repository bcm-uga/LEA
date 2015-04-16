test_extension <- function(name, extension)
{
    # obtain the extension of name
    ext = getExtension(basename(name))

    # if not the correct extension, stop
    if (ext != extension) {
        p = paste("'input_file' format and extension have to be \".", 
            extension, "\" (not \".",ext,"\").", sep="")
        stop(p)
    } 

    return(ext);
}
