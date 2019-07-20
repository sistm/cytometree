#' Data transformation using log10
#'
#' @keywords internal

Log10Transformation <- function(M, num_col, a = 1, b = 1){
  
  M[,num_col] <- log10(a*M[,num_col]+b)
  
  return(M)
}
