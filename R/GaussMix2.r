#' trimodal normal mixture distribution.
#'
#' @keywords internal
#'
#' @importFrom stats dnorm
GaussMix2 <- function(x, mu1, mu2, mu3, 
                      sd1, sd2, sd3, 
                      pi1, pi2, pi3)
{
  comp1 <- pi1*stats::dnorm(x, mu1, sd1)
  comp2 <- pi2*stats::dnorm(x, mu2, sd2)
  comp3 <- pi3*stats::dnorm(x, mu3, sd3)
  return(comp1 + comp2 + comp3)
}