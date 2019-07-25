#' Data transformation using asinh
#'
#' @keywords internal

AsinhTransformation <- function(M, num_col, a = 0, b = 0.2, c = 0){
  
  M[,num_col] <- asinh((a+b*M[,num_col])+c)
  
  return(M)
}
