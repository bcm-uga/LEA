write.geno <- function(R, output.file) 
{
    if(missing(R))
        stop("'R' argument is missing.")
    else if (!(is.matrix(R) || is.data.frame(R) || is.vector(R)))
        stop("'R' argument has to be of type matrix, data.frame or vector.")
    else if (is.vector(R))
        R = matrix(R,ncol=1,nrow=length(R))
    else if (is.data.frame(R))
        R = as.matrix(R)

    output.file = test_character("output.file", output.file, NULL)

    R[which(is.na(R))] = 9
    R[which(is.nan(R))] = 9

    if(any(R != 2 & R != 1 & R != 0 & R != 9))
        stop("'R' matrix can only contains 0, 1, 2 or 9.") 

    write.table(t(R), output.file, col.names=FALSE,row.names=FALSE,sep="");
    return(output.file);
}
