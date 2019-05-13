#' MLE and AIC from mixture of 2 normally distributed observations
#
#' @author Boris Hejblum, Anthony Devaux
#'
#' @param x Numeric vector of observations
#' @param init Numeric vector of parameters (p, mu1, mu2, sigma1, sigma2)
#'
#' @import marqLevAlgParallel
#' 
#' @return List with AIC and parameters
#' 
#' @export
#' 
#' @examples
#' xx <- c(rnorm(100,4,2), rnorm(200, -25, 0.1))
#' aic_2_gauss(xx, init = "kmeans")
#' aic_2_gauss(xx, init = c(0.5, 0, -5, 2, 1))
#'

aic_2_gauss <- function(x, init, maxit = 15, ncore = 1){
  
  browser()
    
  if (class(x)!="numeric" & class(x)!="integer"){
    stop("data vector must be numeric !")
  }
  if (init=="kmeans"){
   
    kmeans_resu <- kmeans(x = x, centers = 2)
    p <- kmeans_resu$size[1]/length(x)
    m1 <- kmeans_resu$centers[1]
    m2 <- kmeans_resu$centers[2]
    s1 <- sd(x[kmeans_resu$cluster==1])
    s2 <- sd(x[kmeans_resu$cluster==2])
    
    init <- c(p,m1,m2,s1,s2)
    
  }else{
    if (class(init)!="numeric"){
      stop("init vector must be numeric !")
    }else if (length(init)!=5){
      stop("init must be of length 5 \n(with the followoing parameters: p, mu1, mu2, sigma1, sigma2)")
    }
  } 
  
  #globals
  n <- length(x)
  init[1] <- logit(init[1])
  
  mloglik_2gauss <- function(b){
    
    p <- expit(b[1])
    mu1 <- b[2]
    mu2 <- b[3]
    s1 <- sqrt(b[4]^2)
    s2 <- sqrt(b[5]^2)
    
    indiv_ll <- sapply(x, function(y){
      log(p * exp(-(y-mu1)^2/(2*s1^2))/(s1 * sqrt(2 * pi)) + (1-p) * exp(-(y-mu2)^2/(2*s2^2))/(s2 * sqrt(2 * pi)))
    })
    return(sum(indiv_ll))
    
  }
  
  deriv_p_llh_2gauss <- function(b){
    
  }
  
  deriv_mu1_llh_2gauss <- function(b){
    
  }
  
  deriv_mu2_llh_2gauss <- function(b){
    
  }
  
  deriv_sigma1_llh_2gauss <- function(b){
    
  }
  
  deriv_sigma2_llh_2gauss <- function(b){
    
  }
  
  gradient <- c(deriv_p_llh_2gauss,
                deriv_mu1_llh_2gauss,
                deriv_mu2_llh_2gauss,
                deriv_sigma1_llh_2gauss,
                deriv_sigma2_llh_2gauss)
  
  resu <- marqLevAlgParallel::marqLevAlg(b = init, fn = mloglik_2gauss, maxiter = maxit,
                                         nproc = ncore, minimize = FALSE)
  
  p_opt <- expit(resu$b[1])
  mu1_opt <- resu$b[2]
  mu2_opt <- resu$b[3]
  s1_opt <- sqrt(resu$b[4]^2)
  s2_opt <- sqrt(resu$b[5]^2)
  AIC_opt <- 2*resu$fn.value + 2*5
  
  indiv_ll1 <- sapply(x, function(y){
    log(p_opt) - (y-mu1_opt)^2/(2*s1_opt^2) - log(s1_opt) - log(2*pi)/2
  })
  indiv_ll2 <- sapply(x, function(y){
    log(1-p_opt) -(y-mu2_opt)^2/(2*s2_opt^2) - log(s2_opt) - log(2*pi)/2
  })
  indiv_clustering <- apply(cbind(indiv_ll1, indiv_ll2), 1, which.max)
  
  return(list(
    "p" = p_opt,
    "mu1" = mu1_opt,
    "mu2" = mu2_opt,
    "s1" = s1_opt,
    "s2" = s2_opt,
    "AIC" = AIC_opt,
    "cluster" = indiv_clustering
  ))
  
}
