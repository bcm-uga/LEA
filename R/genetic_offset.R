genetic.offset <-  function(input, 
                            env, 
                            new.env, 
                            pop.labels,
                            K,
                            pca = FALSE,
                            candidate.loci = NULL) {
  
  
  haploidisation <- function(genotype, pop.labels, env, new.env){
    
    # internal function
    #This function aims at generating a new haploid genotype matrix from a diploid genotype matrix
    # It will also generate new matrices of env and new.env that will be twice the size of 
    # the initial matrices.
    
    # store n and L values
    n <- nrow(genotype)
    L <- ncol(genotype)
    
    # Creation of a new matrix with 2*n individuals and L locus
    haploid.genotype <- matrix(rep(0,(2*n*L)), nrow=2*n, ncol=L)
    # New vector of pop.labels
    haploid.pop.labels <- c()
    # Haploid env matrices
    haploid.env <- c()
    haploid.new.env <- c()
    
    # For all diploid individuals, we will generate two diploid individual
    for (i in seq(1,n)){
      
      current.individual <- genotype[i,]
      
      heterozygote.position <- current.individual == 1
      homozygote.position <- current.individual == 2
      random.vec <- rbinom(n=L, size=1, prob=0.5)
      
      # We will create two individuals following this logic : 
      
      # If at a given locus, diploid ind had 2 version of the derived allele
      # Then at this given locus, the value for both haploid ind will be one
      
      # If at a given locus, diploid ind had 1 version of the derived allele
      # Then we will affect randomly the 1 value to either the first or the second haploid ind
      haploid.ind.1 <- homozygote.position + heterozygote.position*random.vec
      haploid.ind.2 <- homozygote.position + (heterozygote.position - heterozygote.position*random.vec)
      
      haploid.genotype[(2*(i-1) + 1),] <- haploid.ind.1
      haploid.genotype[2*i,] <- haploid.ind.2
      
      current.pop <- pop.labels[i]
      haploid.pop.labels <- c(haploid.pop.labels, rep(current.pop, 2))
      
      current.env <- env[i,]
      haploid.env <- rbind(haploid.env, current.env, current.env)
      
      current.new.env <- new.env[i,]
      haploid.new.env <- rbind(haploid.new.env, current.new.env, current.new.env)
    }
    
    return(list(haploid.genotype=haploid.genotype, haploid.pop.labels=haploid.pop.labels, haploid.env=haploid.env, haploid.new.env=haploid.new.env))
  }
  
            
            ## Check input response matrix 
            ## LEA  
            if (is.character(input)){
              warning("Reading large input files with 'read.lfmm()' may be slow. See 'data.table::fread()' for fast import.")
              Y <- read.lfmm(input)
              lst.unique <- unique(as.numeric(Y))
              if (9 %in% lst.unique){
                stop("'input' file contains missing data (9's). Use the 'impute()' function to impute them.")
              }
              if (-9 %in% lst.unique){
                stop("'input' file contains missing data (-9's). Use the 'impute()' function to impute them.")
              }
            } else {
              ## Y is an R object       
              if (is.null(input)){
                stop("NULL value for argument 'input'.")
              }
              Y <- as.matrix(input)
              Y[Y == 9] <- NA
              Y[Y == -9] <- NA
              lst.unique <- unique(as.numeric(Y))
              if (anyNA(Y)) {
                stop("The input matrix contains missing values: NA, 9 or -9 not allowed. Use the 'write.geno()' and 'impute()' functions to impute them.")
              }
            }
            
            
            ## Check independent/covariate env matrix  
            ## LEA 
            if (is.character(env)){
              X <- read.env(env)
              if (anyNA(X)){
                stop("'env' file contains missing data (NA).")
              }
            } else {
              if (is.null(env)){
                stop("NULL value for argument 'env'.")
              }
              X <- as.matrix(env)
              if (anyNA(X)) {
                stop("The environmental matrix contains NA.")
              }
            }
            
            ## Check new.env matrix  
            ## LEA 
            if (is.character(new.env)){
              X.new <- read.env(new.env)
              if (anyNA(X.new)){
                stop("'new.env' file contains missing data (NA).")
              }
            } else {
              if (is.null(new.env)){
                stop("NULL value for argument 'new.env'.")
              }
              X.new <- as.matrix(new.env)
              if (anyNA(X.new)) {
                stop("The new environmental matrix contains NA.")
              }
            }
            
            ## Checking ploidy
            # cat("Checking ploidy.","\n")
            if (max(lst.unique) > 2) stop("Only haploid or diploid genomes are allowed. For polypoids, perform haploidization (phasing) 
                                          and create copies of individual environmental data.")
            
            ### HAPLOIDIZATION: RANDOM PHASE
            ### Random phase and duplication (n -> 2n)
            
            if (max(lst.unique) == 2){
              cat("Random phase and duplication of samples: 
                  Creating two haploid genomes from each individual genome.","\n")
              haplo <- haploidisation(Y, pop.labels, X, X.new)
              Y <- haplo$haploid.genotype
              pop.labels <- haplo$haploid.pop.labels
              X <- haplo$haploid.env
              X.new <- haplo$haploid.new.env
            }
            
            
            d <-  ncol(X) #number of environmental variables
            d.new <-  ncol(X.new) #number of environmental variables
            if (d.new != d){
              stop("Number of columns in 'new.env' matrix is not equal to the number of columns in 'env' matrix")    
            }
            
            
            
            n <-  nrow(X) #number of individuals
            n.new <-  nrow(X.new) #number of individuals
            if (n.new != n){
              stop("Number of rows in 'new.env' matrix is not equal to the number of rows in 'env' matrix")    
            }
            
            if (nrow(Y) != n){
              stop("Number of rows in the input matrix not equal to the number of rows in the 'env' matrix")    
            }
            
            if (n < d) {
              stop("The environmental covariate matrix contains more columns (d) than rows (n).")
            }
            
            ## Check candidate.loci
            if (is.null(candidate.loci)){
              candidate.loci <- seq(1, ncol(Y))
            } else {
              if (sum(!as.numeric(candidate.loci)) > 0 | max(candidate.loci) > ncol(Y)){
                stop("candidate.loci must be encoded as numeric values not exceeding the total number of columns in the genotype matrix.")}
            }
            
            ## Check pop.labels
            if (length(pop.labels) != n){
              stop("Number of pop labels not equal to the number of rows in the input matrix")    
            }
            
            unique.labels <- sort(unique(pop.labels))
            tab <- sapply(unique.labels, FUN = function(x) sum(pop.labels == x))
            if (min(tab) == 1){
              stop("Population samples in pop.labels must have more than one individual.")    
            }
            
            object <- lfmm2(input = Y, env = X, K = K)
            
            ## compute effect sizes (B matrix)
            U <- object@U
            mod.lm <- lm(Y ~ ., data = data.frame(X, U)) 
            sm <- summary(mod.lm)
            effect.size <- sapply(sm, FUN = function(x) x$coeff[1:(K + d + 1), 1])
            rm(sm)
            
            X.exp <- cbind(rep(1.0, n), X, U)
            Y.fit <- X.exp %*% effect.size
            
            X.pred <- cbind(rep(1.0, n), X.new, U)
            Y.pred <- X.pred %*% effect.size
            
            L <- ncol(Y.fit)
            
            ### Le calcul des offsets s'appuie sur l'ACP  de matrice Y et non sur la matrice Yst (deux populations)
            ### implementer l'option 'coefficient de determination' 
            ##  if (!pca) 'coefficient de determination' else
            
            l.fit <- NULL
            l.pred <- NULL
            l.fitpred <- NULL
            pop.labels.double <- c(pop.labels,pop.labels)
            Y.fit.pred <- rbind(Y.fit, Y.pred)
            offset.r2 <- c()
            
            
            for (i in unique.labels){
              if (!pca){
                genome_current_pop <- Y.fit.pred[pop.labels.double==i,]
                
                #### Preparing data for the regression
                # Count number of individuals
                nb.ind <- dim(genome_current_pop)[1]/2
                # Labels for current and adapted population
                pop.labels.current.futur <- c(rep(1, nb.ind), rep(2, nb.ind))
                genome_current_pop <- scale(genome_current_pop)
                
                regression <- summary(lm(genome_current_pop ~ as.factor(as.numeric(pop.labels.current.futur))))
                fst <- sapply(regression, FUN = function(smr) smr$r.squared)
                offset.r2 <- c(offset.r2, mean(fst[candidate.loci]))
              }else{
                l.fit <- c(l.fit, prcomp(Y.fit[pop.labels == i, candidate.loci], scale = TRUE)$sdev[1]^2)
                l.pred <- c(l.pred, prcomp(Y.pred[pop.labels == i, candidate.loci], scale = TRUE)$sdev[1]^2)
                l.fitpred <- c(l.fitpred, prcomp(rbind(Y.fit[pop.labels == i, candidate.loci], Y.pred[pop.labels == i, candidate.loci]), 
                                                 scale = TRUE)$sdev[1]^2)
              }
            }
            
            if (!pca){
              names(offset.r2) <- unique.labels
              return(offset.r2)
            } else {
              # 1 - F_LT = (1 - F_ST)/(1 - F_SL)
              offset <- 1 - (L - l.fitpred)/(L - (l.pred + l.fit)/2)
              
              offset <- rbind(l.fitpred/L, offset,  l.fit/L, l.pred/L)
              colnames(offset) <- unique.labels
              rownames(offset) <- c("pca.offset", "Flt.offset", "F.fit","Fst.offset")
              
              return(offset)
            }

}



