#'@param 
#'
#'@author Chariff Alkhassim
#'
#'@export 
#'

bootstrapCI <- function(stat,n,alpha){
  x_bar <- mean(stat)
  m <- length(stat)
  x_bar_star <- rowMeans(matrix(sample(stat,m*n,replace=TRUE),ncol=m))
  delta_star <- sort(x_bar_star) - x_bar
  percentile_sup <- floor(n*(1-alpha))
  percentile_inf <- n - percentile_sup
  inf_val <- x_bar - delta_star[percentile_sup]
  sup_val <- x_bar - delta_star[percentile_inf]
  return(list("inf"= inf_val,"sup"= sup_val))
}



