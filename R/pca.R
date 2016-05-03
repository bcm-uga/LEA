pca <- function(input.file, 
    K, 
    center = TRUE, 
    scale  = FALSE) 
{
    # test arguments and init
    # input file
    input.file = test_character("input.file", input.file, NULL)
    # check extension and convert if necessary
    input.file = test_input_file(input.file, "lfmm")
    input.file = normalizePath(input.file)
    #K
    K = test_integer("K", K, 0);
    if (K < 0)
        stop("'K positive.")
    # center
    center = test_logical("center", center, 1);
    # scaled
    scale = test_logical("scale", scale, 0);
    # dir 
    dir = setExtension(paste(dirname(normalizePath(input.file)), "/",
        basename(input.file), sep=""), ".pca/")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    dir = normalizePath(dir)

    tmp = basename(setExtension(basename(input.file), ""))
    # eigenvalues file 
    eigenvalue.file = paste(dir, "/", tmp, ".eigenvalues", sep="")
    # eigenvectors file 
    eigenvector.file = paste(dir, "/", tmp, ".eigenvectors", sep="")
    # standard deviation file 
    sdev.file = paste(dir, "/", tmp, ".sdev", sep="")
    # x file 
    projection.file = paste(dir, "/", tmp, ".projections", sep="")

    print("******************************");
    print(" Principal Component Analysis ");
    print("******************************");

    # run
    L = 0; n = 0;
        resC = .C("R_pca", 
        as.character(input.file),
        as.character(eigenvalue.file),
        as.character(eigenvector.file),
        as.character(sdev.file),
        as.character(projection.file),
        n = as.integer(n),
        L = as.integer(L),
        K = as.integer(K),
        as.integer(center),
        as.integer(scale)
    );

    # save res 
    res = new("pcaProject");
    res@projDir = paste(dirname(normalizePath(input.file)), "/", sep = "")
    res@pcaDir = paste(basename(setExtension(basename(input.file), ".pca/")), 
        "/", sep = "")
    res@n = as.integer(resC$n);
    res@L = as.integer(resC$L);
    res@K = as.integer(resC$K);
    res@center = center;
    res@scale = scale;
    res@input.file = basename(input.file);
    res@eigenvalue.file = basename(eigenvalue.file);
    res@eigenvector.file = basename(eigenvector.file);
    res@sdev.file = basename(sdev.file);
    res@projection.file = basename(projection.file);

    res@pcaProject.file = basename(setExtension(basename(input.file), ".pcaProject"))
    save.pcaProject(res, res@pcaProject.file); 

    res
}
