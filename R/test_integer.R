test_integer <- function(name, param, default)
{
    if(missing(param)) {
        if(is.null(default)) {
            p = paste("'",name,"' argument is missing.", sep="");
            stop(p)
        } else 
            return(default);
    } else {
        if(is.double(param))
            param = as.integer(param)
            
        if(!is.integer(param)) {
            p = paste("'",name,"' argument has to be of type integer.", sep="");
            stop(p);
        }
    }
    return(param);
}
