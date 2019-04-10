#' Gridsearch
#' 
#' @export
#' 
#' @examples
#' x <- c(rnorm(100,4,2), rnorm(200, -25, 0.1))
#' resu_aic_1_gauss <- aic_1_gauss(x)
#' res <- Gridsearch(x, resu_aic_1_gauss)


Gridsearch <- function(x, resu_aic_1_gauss, iter_max = 15, ntry = 30, mixture = 2){

  #browser()
    
  if (mixture == 2){
    list_aic <- lapply(1:ntry, function(i){
      
      init_p <- 0.5
      init_mu1 <- rnorm(1, resu_aic_1_gauss$mu, resu_aic_1_gauss$var_mu)
      init_mu2 <- rnorm(1, resu_aic_1_gauss$mu, resu_aic_1_gauss$var_mu)
      init_sigma1 <- rnorm(1, resu_aic_1_gauss$sigma, resu_aic_1_gauss$var_sigma)
      init_sigma2 <- rnorm(1, resu_aic_1_gauss$mu, resu_aic_1_gauss$var_mu)
      return(aic_2_gauss(x, init = c(init_p, 
                                     init_mu1, init_mu2,
                                     init_sigma1, init_sigma2), maxit = iter_max))
    })
  }else if (mixture == 3){
    list_aic <- lapply(1:ntry, function(i){
      
      init_p1 <- 1/3
      init_p2 <- 1/3
      init_mu1 <- rnorm(1, resu_aic_1_gauss$mu, resu_aic_1_gauss$var_mu)
      init_mu2 <- rnorm(1, resu_aic_1_gauss$mu, resu_aic_1_gauss$var_mu)
      init_mu3 <- rnorm(1, resu_aic_1_gauss$mu, resu_aic_1_gauss$var_mu)
      init_sigma1 <- rnorm(1, resu_aic_1_gauss$sigma, resu_aic_1_gauss$var_sigma)
      init_sigma2 <- rnorm(1, resu_aic_1_gauss$mu, resu_aic_1_gauss$var_mu)
      init_sigma3 <- rnorm(1, resu_aic_1_gauss$mu, resu_aic_1_gauss$var_mu)
      return(aic_3_gauss(x, init = c(init_p1, init_p2, 
                                     init_mu1, init_mu2, init_mu3,
                                     init_sigma1, init_sigma2, init_sigma3), maxit = iter_max))
    })
  }
  
  best_aic <- list_aic[[which.min(sapply(list_aic, "[[", "AIC"))]]
}