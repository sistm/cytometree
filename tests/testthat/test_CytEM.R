library(cytometree)
data(DLBCL)
cellevents <- DLBCL[, c("FL1", "FL2", "FL4")]

context("test CytEM")
myminleaf <- 100
test_that("minleaf works as expected", {
  expect_gte(length(CytEM(M = as.matrix(cellevents), indices = 1:nrow(cellevents), minleaf = myminleaf, level = 1, t = 0.1, force_marker = NULL)[["child"]][["L"]]), myminleaf)
  expect_gte(length(CytEM(M = as.matrix(cellevents), indices = 1:nrow(cellevents), minleaf = myminleaf, level = 1, t = 0.1, force_marker = NULL)[["child"]][["R"]]), myminleaf)
})
