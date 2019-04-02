context("test-aic_1_gauss")

x <- rnorm(10, m=5, sd=2)
res <- aic_1_gauss(x)

test_that("loglikelihood is correct", {
  expect_equal(res$loglikelihood, 
               sum(dnorm(x, m = res$mu, sd = res$sigma, log = TRUE))
  )
})
