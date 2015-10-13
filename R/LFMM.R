lfmm <- function(input.file, 
    environment.file, 
    K,
    project = "continue", 
    d = 0,
    all = FALSE,
    missing.data = FALSE,
    CPU = 1,
    iterations = 10000,
    burnin = 5000,
    seed = -1, 
    repetitions = 1,
    epsilon.noise = 1e-3,
    epsilon.b = 1000,
    random.init = TRUE)
{

    ###########################
    # test arguments and init #
    ###########################

    # input file
    input.file = test_character("input.file", input.file, NULL)
    input.file = test_input_file(input.file, "lfmm")
    input.file = normalizePath(input.file)
    # cov file
    environment.file = test_character("environment.file", 
        environment.file, NULL)
    environment.file = normalizePath(environment.file)
    # check extension
    test_extension(environment.file, "env")
    # K
    for (k in 1:length(K)) {
        K[k] = test_integer("K", K[k], NULL)
        if (K[k] < 0)
            stop("'K' argument has to be not negative.")
    }
    # nd
    if (!missing(d)) {
        for (ndd in 1:length(d)) {
            d[ndd] = test_integer("nd", d[ndd], 1)
            if(d[ndd] <= 0)
                d[ndd] = 1;
        }
    } else {
        v = dim(read.env(environment.file))
        nD = v[2]
        d=1:nD
    }
    # all
    all = test_logical("all",all, FALSE)
    # output.file
    output.file = setExtension(input.file, "")
    # missing.data  
    missing.data = test_logical("missing.data", missing.data, FALSE)
    # CPU    
    CPU = test_integer("CPU", CPU, 1)
    if (CPU <= 0)
                CPU = 1;
    if(Sys.info()['sysname'] == "Windows")
        CPU = 1;
        # iterations
    iterations = test_integer("iterations", iterations, 10000)
    if (iterations <= 0)
                stop("'iterations' argument has to be positive.")
        # burnin
    burnin = test_integer("burnin", burnin, 5000)
    if (burnin <= 0)
        stop("'burnin' argument has to be positive.")
    if (burnin >= iterations) {
        stop("the number of iterations for burnin (burnin) is greater",
            "than the number total of iterations (iterations)")
    }
    # repetitions
    repetitions = test_integer("repetitions", repetitions, 1)
    # epsilon.noise
    epsilon.noise = test_double("epsilon.noise", epsilon.noise, 1e-3)
    if (epsilon.noise < 0)
        epsilon.noise = 1e-3;
    # b.epsilon
    epsilon.b = test_double("epsilon.b", epsilon.b, 1000)
    if (epsilon.b < 0)
        epsilon.b = 1000;
    # random.init 
    random.init = test_logical("random.init", random.init, TRUE)
    # project
    if (missing(project)) 
        project = "continue"
    else if (!(project %in% c("continue", "new", "force"))) 
        stop("A project argument can be 'continue', 'new' or 'force'.");


    ####################
    # call the project #
    ####################

    proj = projectLfmmLoad(input.file, environment.file, project)

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

        for (k in K) {

            if (all) {
                print("*************************************");
                p = paste("* lfmm K =",k," repetition",r," all     *");
                print(p);
                print("*************************************");

                # find the number of the run
                if (length(proj@K) > 0)
                    re = length(which(proj@K == k & proj@all == all)) + 1
                else 
                    re = 1

                # create a directory for the run
                tmp  = setExtension(basename(proj@lfmmProject.file), ".lfmm/")
                dir = paste(proj@directory, "K", k, "/run", re, "/", sep="")
                dir.create(dir, showWarnings = FALSE, recursive = TRUE)

                output.prefix = paste(dir, basename(output.file),"_r", re, 
                    sep="")

                dic = 0
                dev = 0
                L = 0
                n = 0
                D = 0
                resC = .C("R_LFMM", 
                    as.character(currentDir(input.file)),
                    as.character(currentDir(output.prefix)),
                    as.character(currentDir(environment.file)),
                    n = as.integer(n),
                    L = as.integer(L),
                    D = as.integer(D),
                    as.integer(d),
                    as.integer(k),
                    as.integer(iterations),
                    as.integer(burnin),
                    as.integer(CPU),
                    s = as.integer(s),
                    as.integer(missing.data),
                    as.integer(all),
                    dic = as.double(dic),
                    dev = as.double(dev),
                    epsilon.noise,
                    epsilon.b,
                    random.init
                );

                for (nd in 1:nD) { 
                    # creation of the res file
                    res = new("lfmmClass");
                    res@directory = dir;
                    res@zscore.file = normalizePath(paste(output.prefix,
                        "_a",nd, ".",k,".zscore",sep=""));
                    res@lfmmClass.file = paste(output.prefix, 
                        "_a",nd,".",k,".lfmmClass",sep="");
                    res@K = as.integer(k);
                    res@run = as.integer(re);
                    res@d = as.integer(nd);
                    res@Niter = as.integer(iterations);
                    res@burn = as.integer(burnin);
                    res@CPU = as.integer(CPU);
                    res@seed = as.integer(s);
                    res@missing.data = missing.data;
                    res@all = all;
                    res@epsilon.noise = epsilon.noise;
                    res@epsilon.b = epsilon.b;
                    res@random.init = random.init;
                    res@seed = resC$s
                    # res@inflationFactor = inflationFactorEstimation(res)
                    res@deviance = resC$dev;
                    res@DIC = resC$dic;
                    s = resC$s
                    save.lfmmClass(res, res@lfmmClass.file)

                    proj@n = as.integer(resC$n);
                    proj@L = as.integer(resC$L);
                    proj@D = as.integer(resC$D);
                    proj = addRun.lfmmProject(proj, res);
                }
            } else {
                for (nd in d) {
                    print("********************************");
                    p = paste("* K =",k," repetition",r," d =",nd,"  *");
                    print(p);
                    print("********************************");
            
                    # find the number of the run
                    if (length(proj@K) > 0)
                        re = length(which(proj@K == k & proj@d == nd 
                            & proj@all == all)) + 1
                    else 
                        re = 1

                    # create a directory for the run
                    tmp  = setExtension(basename(proj@lfmmProject.file), 
                        ".lfmm/")
                    dir = paste(proj@directory, "K", k, "/run", re, "/", sep="")
                    dir.create(dir, showWarnings = FALSE, recursive = TRUE)

                    output.prefix = paste(dir, basename(output.file),"_r", re, 
                        sep="")
    
                    dic = 0; dev = 0; L = 0;
                    n = 0; D = 0;
                    resC =  .C("R_LFMM", 
                        as.character(currentDir(input.file)),
                        as.character(currentDir(output.prefix)),
                        as.character(currentDir(environment.file)),
                        n = as.integer(n),
                        L = as.integer(L),
                        D = as.integer(D),
                        as.integer(nd),
                        as.integer(k),
                        as.integer(iterations),
                        as.integer(burnin),
                        as.integer(CPU),
                        s = as.integer(s),
                        as.integer(missing.data),
                        as.integer(all),
                        dic = as.double(dic),
                        dev = as.double(dev),
                        epsilon.noise,
                        epsilon.b,
                        random.init
                    );

                    # creation of the res file
                    res = new("lfmmClass");
                    res@zscore.file = normalizePath(paste(output.prefix,
                        "_s",nd,".",k,".zscore",sep=""));
                    res@lfmmClass.file = paste(output.prefix, 
                        "_s",nd,".",k,".lfmmClass",sep="");
                    res@directory = getwd();
                    res@K = as.integer(k);
                    res@run = as.integer(re);
                    res@d = as.integer(nd);
                    res@Niter = as.integer(iterations);
                    res@burn = as.integer(burnin);
                    res@CPU = as.integer(CPU);
                    res@seed = as.integer(s);
                    res@missing.data = missing.data;
                    res@all = all;
                    res@epsilon.noise = epsilon.noise;
                    res@epsilon.b = epsilon.b;
                    res@random.init = random.init;
                    res@seed = resC$s
                    # res@inflationFactor = inflationFactorEstimation(res)
                    res@deviance = resC$dev;
                    res@DIC = resC$dic;
                    s = resC$s
                    save.lfmmClass(res, res@lfmmClass.file);

                    proj@n = as.integer(resC$n);
                    proj@L = as.integer(resC$L);
                    proj@D = as.integer(resC$D);
                    proj = addRun.lfmmProject(proj, res);
                    save.lfmmProject(proj)
                } 
            }
        }    
    }

    return(proj);
} 
