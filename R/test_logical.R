test_logical <- function(name, param, default)
{
    if(missing(param)) {
        if(is.null(default)) {
            p = paste("'",name,"' argument is missing.", sep="");
            stop(p)
        } else 
            return(default);
    } else {
        if(!is.logical(param)) {
            p = paste("'",name,"' argument has to be of type logical.", sep="");
            stop(p);
        }
    }
    return(param);
}
