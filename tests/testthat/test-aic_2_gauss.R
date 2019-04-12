context("test-aic_2_gauss")

x <- c(rnorm(100,m=4,sd=2), rnorm(200,m=-25,sd=0.1))
res <- aic_2_gauss(x, init = "kmeans", maxit = 100)

test_that("loglikelihood is correct", {
  expect_equal(res$AIC, 
               -2*sum(
                 log(res$p * exp(-(x-res$mu1)^2/(2*res$s1^2))/(res$s1 * sqrt(2 * pi)) + 
                            (1-res$p) * exp(-(x-res$mu2)^2/(2*res$s2^2))/(res$s2 * sqrt(2 * pi)))
                 ) + 2*5)
})

