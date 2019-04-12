context("test-logit")

x <- 0.15

test_that("result is correct", {
  expect_equal(logit(x), log(x/(1 - x)))
})
