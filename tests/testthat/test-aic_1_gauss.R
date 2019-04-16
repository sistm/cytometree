context("test-aic_1_gauss")

x <- rnorm(1000, m=5, sd=2)
res <- aic_1_gauss(x)

test_that("loglikelihood is correct", {
  expect_equal(res$loglikelihood, 
               sum(dnorm(x, m = res$mu, sd = res$sigma, log = TRUE))
  )
})

test_that("approximation is correct", {
  expect_equal(res$loglikelihood,
               sum(dnorm(x, m = 5, sd = 2, log = TRUE)),
               tolerance = 0.01
  )
})