#' MLE and AIC from normally distributed observations
#
#' @author Anthony Devaux, Boris Hejblum
#'
#' @param x a numeric vector of observations
#' 
#' @return a list with AIC and parameters
#' 
#' @export
#' 
#' @examples
#' x <- rnorm(10, m=5, sd=2)
#' res <- aic_1_gauss(x)
#'
#' res$mu
#' res$sigma
#' res$loglikelihood
#' sum(dnorm(x, m = res$mu, sd = res$sigma, log = TRUE))

aic_1_gauss <- function(x){
  
  if (class(x)!="numeric" & class(x)!="integer"){
    stop("x vector must be numeric !")
  }
  
  n <- length(x)
  
  mu_mle <- mean(x)
  sigma2_mle <- mean((x - mu_mle)^2)
  
  var_mu_FisherInfo <- sigma2_mle/n
  var_sigma_FisherInfo <- (2 * sigma2_mle^2)/n
  
  log_likelihood <- - sum((x-mu_mle)^2) / (2*sigma2_mle) - n/2*log(2*pi*sigma2_mle)
  
  aic <- -2*log_likelihood + 2*2
  
  return(list("mu" = mu_mle, 
              "sigma" = sqrt(sigma2_mle), 
              "var_mu" = var_mu_FisherInfo, 
              "var_sigma" = var_sigma_FisherInfo, 
              "loglikelihood" = log_likelihood,
              "AIC" = aic))
  
}
