cross.entropy.estimation <- function(input.file, K, masked.file, Q.file, 
    G.file, ploidy = 2) {

    # test arguments and init
    # input file
    input.file = test_character("input.file", input.file, NULL)
    # check extension and convert if necessary
    input.file = test_input_file(input.file, "geno")
    # K
    K = test_integer("K", K, NULL)
    if (K <= 0)
        stop("'K' argument has to be positive.")
    # masked data file
    tmp = setExtension(input.file, "_I.geno");
    masked.file = test_character("masked.file", masked.file, tmp)
    # check extension
    test_extension(masked.file, "geno")
    # Q file    
    tmp = paste(setExtension(input.file, ""), "_I.", K, ".Q", sep="")
    Q.file = test_character("Q.file", Q.file, tmp) 
    # check extension
    test_extension(Q.file, "Q")
    # G file
    tmp = paste(setExtension(input.file, ""), "_I.", K, ".G", sep="")
    G.file = test_character("G.file", G.file, tmp) 
    # check extension
    test_extension(G.file, "G")
    # ploidy
    ploidy = test_integer("ploidy", ploidy, 2)
    if (ploidy <= 0)
        stop("'ploidy' argument has to be positive.")

    print("*************************************");
    print("*    cross-entropy estimation       *");
    print("*************************************");

    # run method
    all.ce = 0;
    masked.ce = 0;
        res = .C("R_crossEntropy", 
        as.character(input.file),
        as.character(masked.file),
        as.character(Q.file),
        as.character(G.file),
        as.integer(K),
        as.integer(ploidy),
        all.ce = as.double(all.ce),
        masked.ce = as.double(masked.ce)
    )

    return(list(masked.ce=res$masked.ce, all.ce=res$all.ce))
}
