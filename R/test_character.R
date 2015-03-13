test_character <- function(name, param, default)
{
    if(missing(param)) {
        if(is.null(default)) {
            p = paste("'",name,"' argument is missing.", sep="");
            stop(p)
        } else 
            return(default);
    } else {
        if(!is.character(param)) {
            p = paste("'",name,"' argument has to be of type character.", 
                sep="");
            stop(p);
        }
    }
    return(param);
}
