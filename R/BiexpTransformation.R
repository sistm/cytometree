#' Data transformation using biexp
#'
#' @importFrom GoFKernel inverse
#'
#' @keywords internal

BiexpTransformation <- function(M, num_col, a = 0.5, b = 1, c = 0.5, d = 1, f = 0, w = 0){
  
  BiexpFunction <- function(x){
    a*exp(b*(x-w))-c*exp(-d*(x-w))+f
  }
  
  InverseFunction <- GoFKernel::inverse(BiexpFunction)
  
  M[,num_col] <- sapply(M[,num_col], InverseFunction)
  
  return(M)
}