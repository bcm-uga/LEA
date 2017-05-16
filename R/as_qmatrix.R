as.qmatrix <- function(Q){
  if (class(Q) != "matrix") Q <- as.matrix(Q)
  if (min(Q) < 0) stop("Q contains negative elements.")
  sumofq <- apply(Q, MARGIN = 1, sum)
  sumofq <- round(sum(sumofq))
  if (sumofq != nrow(Q)) stop("Input matrix is not an ancestry matrix: The sum of ancestry coefficients is not equal to one")
  class(Q) = c("Qmatrix","matrix") 
  return(Q)
}
