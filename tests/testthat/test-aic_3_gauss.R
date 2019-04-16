context("test-aic_3_gauss")

x <- c(rnorm(100,4,2), rnorm(600, -25, 0.1), rnorm(300, 30, 4))
res <- aic_3_gauss(x, init = "kmeans", maxit = 100)
#res <- aic_3_gauss(x, init = c(0.3, 0.3, 0, -10, 10, 2, 2, 2), maxit = 100)

test_that("AIC is correct", {
  expect_equal(res$AIC,
               -2*sum(
                 log(res$p1 * exp(-(x-res$mu1)^2/(2*res$s1^2))/(res$s1 * sqrt(2 * pi)) + 
                       res$p2 * exp(-(x-res$mu2)^2/(2*res$s2^2))/(res$s2 * sqrt(2 * pi))  +
                       res$p3 * exp(-(x-res$mu3)^2/(2*res$s3^2))/(res$s3 * sqrt(2 * pi)))
               ) + 2*8)
})

test_that("approximation is correct", {
  expect_equal(res$AIC,
               -2*sum(log(1/10 * exp(-(x-4)^2/(2*2^2))/(2 * sqrt(2 * pi)) + 
                            6/10 * exp(-(x--25)^2/(2*0.1^2))/(0.1 * sqrt(2 * pi))  +
                            3/10 * exp(-(x-30)^2/(2*4^2))/(4 * sqrt(2 * pi)))) + 2*8,
               tolerance = 0.01
  )
})