#' Compute parameters and AIC from observations with normal distribution
#
#' @author Anthony DEVAUX
#'
#' @param x Numeric vector of observations
#'
#' @importFrom stats 
#' 
#' @return List with AIC and parameters
#' 
#'@examples
#'data_obs <- rnorm(100,2,4)
#'
#'res <- aic_1_gauss(data_obs)
#'
#'res$mu
#'res$aic

aic_1_gauss <- function(x){
  if (class(x)!="numeric"){
    stop("x vector can be numeric !")
  }
  
  n <- length(x)
  
  mu <- sum(x)/n
  sigma <- sqrt(sum((x-mu)^2)/n)
  
  var_mu <- (sigma^2)/n
  var_sigma <- (2*sigma^4)/n
  
  vraisemblance <- prod(((sigma*sqrt(2*pi))^-1)*exp((-(x-mu)^2)/sigma^2))
  
  aic <- -2*log(vraisemblance)+4
  
  theta <- list(mu = mu, sigma = sigma, var_mu = var_mu, var_sigma = var_sigma, aic = aic)
  
  return(theta)
  
}