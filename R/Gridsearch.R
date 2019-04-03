#' Gridsearch
#' 
#' @export
#' 
#' @example 
#' x <- c(rnorm(100,4,2), rnorm(200, -25, 0.1))
#' resu_aic_1_gauss <- aic_1_gauss(x)
#' res <- Gridsearch(x, resu_aic_1_gauss)


Gridsearch <- function(x, resu_aic_1_gauss, iter_max = 15, ntry = 30){
  
  list_aic2 <- lapply(1:ntry, function(i){
    
    init_p <- 0.5
    init_mu1 <- rnorm(1, resu_aic_1_gauss$mu, resu_aic_1_gauss$var_mu)
    init_mu2 <- rnorm(1, resu_aic_1_gauss$mu, resu_aic_1_gauss$var_mu)
    init_sigma1 <- rnorm(1, resu_aic_1_gauss$sigma, resu_aic_1_gauss$var_sigma)
    init_sigma2 <- rnorm(1, resu_aic_1_gauss$mu, resu_aic_1_gauss$var_mu)
    return(aic_2_gauss(x, init = c(init_p, 
                                   init_mu1, init_mu2,
                                   init_sigma1, init_sigma2), maxit = iter_max))
  })
  best_aic2 <- list_aic2[[which.min(sapply(list_aic2, "[[", "AIC"))]]
}