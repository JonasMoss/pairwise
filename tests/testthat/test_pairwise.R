source("aaa.R")

obj <- pairwise(data)
testthat::expect_gt(obj$sigma, 0.89)
testthat::expect_equal(coef.pairwise(obj, reference = NULL), obj$beta)
