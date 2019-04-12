context("test-expit")

x <- -1.8

test_that("result is correct", {
  expect_equal(expit(x), 1/(1 + exp(-x)))
})
