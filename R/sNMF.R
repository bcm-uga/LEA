snmf <- function(input.file, 
    K, 
    project = "continue",
    repetitions = 1,
    CPU = 1, 
    alpha = 10, 
    tolerance = 0.00001, 
    entropy = FALSE,
    percentage = 0.05,
    I, 
    iterations = 200, 
    ploidy = 2, 
    seed = -1, 
    Q.input.file)
{

    ###########################
    # test arguments and init #
    ###########################

    # input file
    input.file = test_character("input.file", input.file, NULL)
    # check extension and convert if necessary 
    input.file = test_input_file(input.file, "geno")
    input.file = normalizePath(input.file)
    # K
    for (k in 1:length(K)) {
        K[k] = test_integer("K", K[k], NULL)
        if (K[k] <= 0)
            stop("'K' argument has to be positive.")
    }
    # alpha
    alpha = test_double("alpha", alpha, 10)
    if (alpha < 0)
        alpha = 0
    # tolerance
    tolerance = test_double("tolerance", tolerance, 0.0001)
    if (tolerance <= 0)
        tolerance = 0.0001
    # entropy
    entropy = test_logical("entropy", entropy, FALSE)
    # percentage
    percentage = test_double("percentage", percentage, 0)
    if (entropy && (percentage < 0 || percentage >= 1))
        percentage = 0.05
    else if (!entropy)
        percentage = 0
    # iterations
    iterations = test_integer("iterations", iterations, 200)
    if (iterations <= 0)
        iterations = 200;
    # ploidy
    ploidy = test_integer("ploidy", ploidy, 0)
    if (ploidy <= 0)
        ploidy = 0;
    # CPU    
    CPU = test_integer("CPU", CPU, 1)
    if (CPU <= 0)
        CPU = 1;
    if(Sys.info()['sysname'] == "Windows")
        CPU = 1;
    # input Q
    Q.input.file = test_character("Q.input.file", Q.input.file, "")
    # test extension
    if (Q.input.file != "")
    test_extension(Q.input.file, "Q")
    # I    
    I = test_integer("I", I, 0)
    if (I < 0)
        stop("'I' argument has to be of type positive.")
    # repetitions
    repetitions = test_integer("repetitions", repetitions, 1)
    # project
    if (missing(project))
    project = "continue"
    else if (!(project %in% c("continue", "new", "force")))
        stop("A project argument can be 'continue', 'new' or 'force'.");
    
    ####################
    # call the project #
    ####################

    proj = projectSnmfLoad(input.file, project) 
    
    ################################
    # launch each run sequentially #
    ################################

    for (r in 1:repetitions) {
        # set the seed
        if (is.na(seed[r]))
            s = -1
        else
            s = seed[r]
        s = test_integer("seed", s, as.integer(runif(1)*.Machine$integer.max))
        if (s == -1)   
            s = as.integer(runif(1)*.Machine$integer.max)
        set.seed(s) # init seed

        # create.dataset 
        if (entropy) {
            masked.file = setExtension((paste(proj@directory, "masked/", 
            basename(input.file), sep="")), "_I.geno")
            masked.file = create.dataset(input.file, masked.file, s, 
                percentage); 
        } else {
            masked.file = input.file
        }
        for (k in K) {
            print("*************************************");
            p = paste("* sNMF K =",k," repetition",r,"     *");
            print(p);
            print("*************************************");

            re = length(which(proj@K == k)) + 1

            # create a directory for the run
            tmp  = basename(setExtension(basename(input.file), ""))
            dir = paste(proj@directory, "K", k, "/run", re, "/", sep="")
            dir.create(dir, showWarnings = FALSE, recursive = TRUE)        

            # Q file
            Q.output.file = paste(dir, tmp, "_r", re ,".",k, ".Q", sep="")
            # G file
            G.output.file = paste(dir, tmp, "_r", re ,".",k, ".G", sep="")

            # TODO on peut aussi tester que le fichier n est pas déjà 
            # existant 
            snmfClass.file = paste(dir, tmp, "_r", re ,".",k, ".snmfClass", 
                sep="")

            all.ce = 0;
            masked.ce = 0;
            n = 0;
            L = 0;
            resC = .C("R_sNMF", 
                as.character(masked.file),
                as.integer(k),
                as.double(alpha),
                as.double(tolerance),
                as.double(0.0),
                as.integer(iterations),
                s = as.integer(s),
                as.integer(ploidy),
                as.integer(CPU),
                as.character(Q.input.file),
                as.character(Q.output.file),
                as.character(G.output.file),
                as.integer(I),
                all.ce = as.double(all.ce),
                masked.ce = as.double(masked.ce),
                n = as.integer(n),
                L = as.integer(L)
            );
    
            # calculate crossEntropy
            if (entropy) {
                ce = cross.entropy.estimation(input.file, k, masked.file,  
                    Q.output.file, G.output.file, ploidy)
                all.ce = ce$all.ce
                masked.ce = ce$masked.ce
            }
        
            # creation of the res file
            res = new("snmfClass")
            res@directory = dir

            # file snmfClass
            res@snmfClass.file  = snmfClass.file;
            res@K = as.integer(k);
            res@run = as.integer(re);
            res@CPU = as.integer(CPU);
            res@seed = resC$s;
            res@alpha = alpha;
            res@percentage = percentage;
            res@I = as.integer(I);
            res@iterations = as.integer(iterations);
            res@entropy = entropy;
            res@tolerance = tolerance;
            res@crossEntropy = masked.ce;
            res@ploidy = as.integer(ploidy);
            res@Q.input.file = Q.input.file;
            res@Q.output.file = normalizePath(Q.output.file);
            res@G.output.file = normalizePath(G.output.file);
            save.snmfClass(res, res@snmfClass.file)

            proj@n = resC$n;
            proj@L = resC$L;
            proj = addRun.snmfProject(proj, res);
            save.snmfProject(proj)
        }
    }

    return(proj);
}
