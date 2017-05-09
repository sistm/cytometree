GaussMix <- function(x, mu1, mu2, sd1, sd2, pi1, pi2)
{
  comp1 <- pi1*dnorm(x, mu1, sd1)
  comp2 <- pi2*dnorm(x, mu2, sd2)
  return(comp1 + comp2)
}
