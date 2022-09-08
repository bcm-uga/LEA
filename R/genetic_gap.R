genetic.gap <-  function(input, 
                         env,
                         new.env = NULL, 
                         pred.env,
                         K = NULL,
                         scale = FALSE,
                         candidate.loci = NULL) {
  
  ## Check input response matrix 
  ## LEA  
  if (is.character(input)){
    warning("Loading large input files with 'read.lfmm()' may be slow. See 'data.table::fread()' for fast import.")
    Y <- read.lfmm(input)
    lst.unique <- unique(as.numeric(Y))
    if (9 %in% lst.unique){
      stop("'input' file contains missing data (9's). Use the 'impute()' function to impute them.")
    }
    if (-9 %in% lst.unique){
      stop("'input' file contains missing data (-9's). Use the 'impute()' function to impute them.")
    }
  } else {
    ## input is an R object       
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
      stop("'env' file contains missing data (NA's).")
    }
  } else {
    if (is.null(env)){
      stop("NULL value for argument 'env'.")
    }
    X <- as.matrix(env)
    if (anyNA(X)) {
      stop("The environmental matrix contains NA's.")
    }
  }
  
  
  ## Check new.env matrix  
  ## LEA 
  if (is.null(new.env)){
    X.new <- X
  } else {
  X.new <- as.matrix(new.env)
  if (anyNA(X.new)) {
    stop("The new environmental matrix contains NA's.")}
  }
  
  ## Check pred.env matrix  
  ## LEA 
  if (is.character(pred.env)){
    X.pred <- read.env(pred.env)
    if (anyNA(X.pred)){
      stop("The 'pred.env' environmental matrix file contains missing data (NA's).")
    }
  } else {
    if (is.null(pred.env)){
      stop("NULL value for argument 'pred.env'.")
    }
    X.pred <- as.matrix(pred.env)
    if (anyNA(X.pred)) {
      stop("The predicted environmental matrix contains NA's.")
    }
  }
  
  
  
  d <- ncol(X) #number of environmental variables
  d.new <-  ncol(X.new) #number of environmental variables
  
  if (d.new != d){
    stop("Number of columns in 'new.env' matrix is not equal to the number of columns in 'env' matrix")    
  }
  d.pred <-  ncol(X.pred) #number of environmental variables
  if (d.pred != d){
    stop("Number of columns in 'pred.env' matrix is not equal to the number of columns in 'env' matrix")    
  }
  
  n <-  nrow(X) #number of individuals
  if (nrow(Y) != n){
    stop("Number of rows in the input matrix not equal to the number of rows in the 'env' matrix")    
  }
  
  if (n < d) {
    stop("The environmental covariate matrix contains more columns (d) than rows (n).")
  }
  
  
  n.new <-  nrow(X.new) #number of test points (individuals) in the new environmental matrix
  n.pred <-  nrow(X.pred) #number of test points (individuals) in the predicted environmental matrix
  
  if (n.new != n.pred){
    stop("Number of rows in 'new.env' matrix is not equal to the number of rows in 'pred.env' matrix")    
  }

  
  ## Check K
  if (is.null(K)){
    stop("Null value for the number of factor in the LFMM.")  
  } 
  
  ##scale option
  if (scale == TRUE){
     m.x <- apply(X, 2, mean)
     sd.x <- apply(X, 2, sd)
     if (sum(sd.x == 0) > 0){
       stop("Error with scale = TRUE: Impossible to standardize a null column (locus).")  
     } 
     X <- t(t(X) - m.x) %*% diag(1/sd.x)
     X.new <- t(t(X.new) - m.x) %*% diag(1/sd.x)
     X.pred <- t(t(X.pred) - m.x) %*% diag(1/sd.x)
  }   
  
  ## Check candidate.loci
  if (is.null(candidate.loci)){
    candidate.loci <- seq(1, ncol(Y))
  } else {
    if (sum(!as.numeric(candidate.loci)) > 0 | max(candidate.loci) > ncol(Y)){
      stop("candidate.loci must be encoded as numeric values not exceeding the total number of columns in the genotype matrix.")}
  }
  
  object <- lfmm2(input = Y, env = X, K = K)
  
  ## compute effect sizes (B matrix)
  U <- object@U
  mod.lm <- lm(Y ~ ., data = data.frame(X, U)) 
  sm <- summary(mod.lm)
  B <- sapply(sm, FUN = function(x) x$coeff[2:(d + 1), 1])
  rm(sm)
  
  M = (X.new - X.pred)  %*% B
  D = diag(M %*% t(M))/ncol(B) 
  
  eig <- eigen(cov(t(B)))
  
  return(list(offset = D, distance = sqrt(D), eigenvalues = eig$values, vectors = eig$vectors))
  
}



