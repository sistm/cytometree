library(cytometree)
data(DLBCL)
cellevents <- DLBCL[, c("FL1", "FL2", "FL4")]

context("test CytEM")

test_that("minleaf works as expected", {
  expect_gte(length(CytEM(M = as.matrix(cellevents), indices = 1:nrow(cellevents), minleaf = 100, level = 1, t = 0.1, force_marker = NULL)[["child"]][["L"]]), minleaf)
  expect_gte(length(CytEM(M = as.matrix(cellevents), indices = 1:nrow(cellevents), minleaf = 100, level = 1, t = 0.1, force_marker = NULL)[["child"]][["R"]]), minleaf)
})
