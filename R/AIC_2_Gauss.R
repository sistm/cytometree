#' Compute parameters and AIC from observations with mixture of 2 normal distribution
#
#' @author Boris Hejblum, Anthony DEVAUX
#'
#' @param x Numeric vector of observations
#' @param init Numeric vector of parameters (p, mu1, mu2, sigma1, sigma2)
#'
#' @importFrom stats marqLevAlg
#' 
#' @return List with AIC and parameters
#' 
#'@examples
#'xx <- c(rnorm(100,4,2), rnorm(200, -25, 0.1))
#'aic_2_gauss(xx, init = c(0.5, 0, 0, 2, 2))

aic_2_gauss <- function(x, init){
  
  if (class(x)!="numeric"){
    stop("data vector can be numeric !")
  }
  
  #globals
  n <- length(x)
  init[1] <- init[1]/(1 - init[1])
  
  mlogvrais_2gauss <- function(b){
    
    p <- 1/(1 + exp(b[1]))
    mu1 <- b[2]
    mu2 <- b[3]
    s1 <- sqrt(b[4]^2)
    s2 <- sqrt(b[5]^2)
    
    indiv_lv <- sapply(x, function(y){
      log(p * exp(-(y-mu1)^2/s1^2)/(s1 * sqrt(2 * pi)) + (1-p) * exp(-(y-mu2)^2/s2^2)/(s2 * sqrt(2 * pi)))
    })
    return(-sum(indiv_lv))
    
  }
  
  resu <- marqLevAlg(b = init, fn = mlogvrais_2gauss)
  return(list(
    "p" = 1/(1 + exp(resu$b[1])),
    "mu1" = resu$b[2],
    "mu2" = resu$b[3],
    "s1" = sqrt(resu$b[4]^2),
    "s2" = sqrt(resu$b[5]^2),
    "AIC" = 2*resu$fn.value + 2*5
  ))
  
}