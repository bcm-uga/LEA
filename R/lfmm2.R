# create class lfmm2
setClass("lfmm2Class",
         slots = c(K = "integer", 
                   lambda = "numeric",
                   U = "matrix",
                   V = "matrix"
                   )
)

lfmm2 <- function(input,
                  env, 
                  K, 
                  lambda = 1e-5){
  
## Check response input matrix 
## LEA  
    if (is.character(input)){
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
      if (anyNA(Y)) {
        stop("The input matrix contains missing values: NA, 9 or -9 not allowed.")
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

  if (length(K) > 1){
    stop("Multiple values of K not allowed.")
  }
  if (lambda <= 0){
    stop("The ridge regularization parameter must be positive.")
  }
  
  d <-  ncol(X) #number of environmental variables
  n <-  nrow(X) #number of individuals
  
  if (nrow(Y) != n){
    stop("Number of rows in the input matrix not equal to the number of rows in the 'env' matrix")    
  }
  
  if (n < d) {
    stop("The environmental covariate matrix X contains more columns (d) than rows (n).")
  }

# run SVD of X: X = Q Sigma R
  
    svx <- svd(x = scale(X, scale = FALSE), nu = n)
    Q <- svx$u
    
    d_lambda <- c(sqrt(lambda/(lambda + svx$d)), rep(1, n-d))
    d_lambda_inv <- c(sqrt((lambda + svx$d)/lambda), rep(1, n-d))
    D_inv <- diag(d_lambda_inv)
    D  <- diag(d_lambda)
 
# run SVD of modified Y    
    svk <- svd(D %*% t(Q) %*% scale(Y, scale = FALSE), nu = K)

# compute the latent matrix W
    W <- Q %*% D_inv %*% tcrossprod(svk$u %*% diag(svk$d[1:K]), svk$v[,1:K])

# compute LFMM factors U and loadings V
# Non orthogonal factors
    U <- crossprod(t(Q %*% D_inv), svk$u %*% diag(svk$d[1:K]))
    #U <- Q %*% D_inv %*% svk$u %*% diag(svk$d[1:K])
    V <- svk$v[,1:K]

    obj <- new("lfmm2Class")
    obj@K <- as.integer(K)
    obj@lambda <- as.numeric(lambda)
    obj@U <- as.matrix(U)
    obj@V <- as.matrix(V)

## LEA 
    return(obj)
}


setGeneric("lfmm2.test", function(object, input, env, 
                                  genomic.control = TRUE, 
                                  linear = TRUE, 
                                  family  = binomial(link = "logit")) matrix);
setMethod("lfmm2.test", "lfmm2Class",
          function(object, 
                   input,
                   env,
                   genomic.control, 
                   linear,
                   family
                   ) {

            ## Check input matrix   
            ## LEA  
            if (is.character(input)){
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
              if (anyNA(Y)) {
                stop("The input matrix contains missing values (NA or 9).")
              }
            }
            
            ## Check independent/covariate matrix  
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
            
            d <-  ncol(X) #number of environmental variables
            n <-  nrow(X) #number of individuals
            
            if (nrow(Y) != n){
              stop("Number of rows in the input matrix not equal to the number of rows in the 'env' matrix")    
            }
            
            if (n < d) {
              stop("The environmental covariate matrix X contains more columns (d) than rows (n).")
            }
   
            p <- ncol(Y)
            p_value <- NULL
            z_score <- NULL
            
            if (linear){
              mod_lm <- lm(Y ~ ., data = data.frame(X, object@U)) 
              sm <- summary(mod_lm)
              p_value <- sapply(sm, FUN = function(x) x$coeff[2:(d + 1), 4])
              z_score <- as.matrix(sapply(sm, FUN = function(x) x$coeff[2:(d + 1), 3]))
            } else {
              for (j in 1:p) {
              mod_glm <- glm(Y[, j] ~ ., data = data.frame(X, object@U), family = family)
              sm <- summary(mod_glm)
              p_value <- rbind(p_value, sm$coeff[2:(d + 1), 4])
              z_score <- rbind(z_score, sm$coeff[2:(d + 1), 3])
              }
            }
            if (genomic.control){
                gif <- apply(z_score^2, 2, median)/qchisq(0.5, df = 1, lower.tail = FALSE)
                p_value <- pchisq(z_score^2/gif, df = 1, lower.tail = FALSE)
            } else {
                gif <- NULL
            }
       res <- list(pvalues = p_value, zscores = z_score, gif = gif)
       return(res)
      }
)


