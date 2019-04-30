#' #' MLE and AIC from mixture of 3 normally distributed observations
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
#' x <- c(rnorm(100,4,2), rnorm(600, -25, 0.1), rnorm(300, 30, 4))
#' res <- aic_3_gauss(x, init = "kmeans", maxit = 100)
#' res <- aic_3_gauss(x, init = c(0.3, 0.3, 0, -10, 10, 2, 2, 2), maxit = 100)
#'

aic_3_gauss <- function(x, init, maxit = 15, ncore = 1){
  
  if (class(x)!="numeric" & class(x)!="integer"){
    stop("data vector must be numeric !")
  }
  if (init == "kmeans"){
    
    kmeans_resu <- kmeans(x = x, centers = 3)
    p1 <- kmeans_resu$size[1]/length(x)
    p2 <- kmeans_resu$size[2]/length(x)
    m1 <- kmeans_resu$centers[1]
    m2 <- kmeans_resu$centers[2]
    m3 <- kmeans_resu$centers[3]
    s1 <- sd(x[kmeans_resu$cluster==1])
    s2 <- sd(x[kmeans_resu$cluster==2])
    s3 <- sd(x[kmeans_resu$cluster==3])
    
    init <- c(p1,p2,m1,m2,m3,s1,s2,s3)
    
  }else{
    if (class(init)!="numeric"){
      stop("init vector must be numeric !")
    }else if (length(init)!=8){
      stop("init must be of length 8 \n(with the followoing parameters: p1, p2, mu1, mu2, m3, sigma1, sigma2, sigma3)")
    }
    if(init[1]>1 | init[1]<0 | init[2]>1 | init[2]<0 | (init[1] + init[2])>1 | (init[1] + init[2])<0){
      stop(" first 2 parameters arte probabilities that must be in [0;1] and sum under 1")
    }
  }
  
  #globals
  n <- length(x)
  init[1] <- logit(init[1])
  init[2] <- logit(init[2])
  
  mloglik_3gauss <- function(b){
    
    p1 <- expit(b[1])
    p2 <- expit(b[2])*(1-p1)
    p3 <- 1- p1 - p2
    mu1 <- b[3]
    mu2 <- b[4]
    mu3 <- b[5]
    s1 <- sqrt(b[6]^2)
    s2 <- sqrt(b[7]^2)
    s3 <- sqrt(b[8]^2)
    
    indiv_ll <- sapply(x, function(y){
      log(p1 * exp(-(y-mu1)^2/(2*s1^2))/(s1 * sqrt(2 * pi)) + 
          p2 * exp(-(y-mu2)^2/(2*s2^2))/(s2 * sqrt(2 * pi))  +
          (1 - p1 - p2) * exp(-(y-mu3)^2/(2*s3^2))/(s3 * sqrt(2 * pi))
      )
    })
    return(sum(indiv_ll))
    
  }
  
  resu <- marqLevAlgParallel::marqLevAlg(b = init, fn = mloglik_3gauss, maxiter = maxit,
                                         nproc = ncore, minimize = FALSE)
  
  p1_opt <- expit(resu$b[1])
  p2_opt <- expit(resu$b[2])*(1 - p1_opt)
  p3_opt <- 1 - p1_opt - p2_opt
  mu1_opt <- resu$b[3]
  mu2_opt <- resu$b[4]
  mu3_opt <- resu$b[5]
  s1_opt <- sqrt(resu$b[6]^2)
  s2_opt <- sqrt(resu$b[7]^2)
  s3_opt <- sqrt(resu$b[8]^2)
  AIC_opt <- 2*resu$fn.value + 2*8
  
  indiv_ll1 <- sapply(x, function(y){
    log(p1_opt) - (y-mu1_opt)^2/(2*s1_opt^2) - log(s1_opt) - log(2*pi)/2
  })
  indiv_ll2 <- sapply(x, function(y){
    log(p2_opt) -(y-mu2_opt)^2/(2*s2_opt^2) - log(s2_opt) - log(2*pi)/2
  })
  indiv_ll3 <- sapply(x, function(y){
    log(1-(p1_opt+p2_opt)) -(y-mu3_opt)^2/(2*s3_opt^2) - log(s3_opt) - log(2*pi)/2
  })
  indiv_clustering <- apply(cbind(indiv_ll1, indiv_ll2, indiv_ll3), 1, which.max)
  
  return(list(
    "p1" = p1_opt,
    "p2" = p2_opt,
    "p3" = p3_opt,
    "mu1" = mu1_opt,
    "mu2" = mu2_opt,
    "mu3" = mu3_opt,
    "s1" = s1_opt,
    "s2" = s2_opt,
    "s3" = s3_opt,
    "AIC" = AIC_opt,
    "cluster" = indiv_clustering
  ))
  
}
