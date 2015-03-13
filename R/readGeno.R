read.geno <- function(input.file) {

    # test arguments and init
    # input file
    if(missing(input.file))
        stop("'input.file' argument is missing.")
    else if (!is.character(input.file))
        stop("'input.file' argument has to be of type character.")
    # check extension 
    test_extension(input.file, "geno")

    x = scan(file = input.file, what = "character", skip = 0, sep ="")

    if(length(x) > 0) {
        M = length(x)
    } else {
        stop("'input.file' is empty.")
    }

    line = strsplit(x[1],NULL)
    N = length(line[[1]])

    return(apply(as.matrix(x),1, function(x){
        as.integer((strsplit(x,NULL))[[1]])}))
}
