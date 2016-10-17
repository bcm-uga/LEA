read.zscore <- function(input.file) {

    # test arguments
    test_character("input.file", input.file, NULL)
    # check extension 
    test_extension(input.file, "zscore")

    R = as.matrix(read.table(input.file));

    return(list(zscores=R[,1],mlog10pvalues=R[,2],pvalues=R[,3]))
}
