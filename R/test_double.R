test_double <- function(name, param, default)
{
    if(missing(param)) {
        if(is.null(default)) {
            p = paste("'",name,"' argument is missing.", sep="");
            stop(p)
        } else 
            return(default);
    } else {
        if(is.integer(param))
            param = as.double(param)
            
        if(!is.double(param)) {
            p = paste("'",name,"' argument has to be of type double.", sep="");
            stop(p);
        }
    }
    return(param);
}
