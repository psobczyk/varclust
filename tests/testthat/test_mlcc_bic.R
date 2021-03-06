context("Testing mlcc.bic")

library(varclust)

test_that("mlcc.bic for small matrix", {
  set.seed(10)
  load("test_data/small_matrix.rda")
  true_segmentation <- rep(1:2, each = 50)
  result <- mlcc.bic(X, numb.clusters = 1:5, max.dim = 2, numb.cores = 1)
  expect_equal(2, result$nClusters)
  scores <- integration(result$segmentation, true_segmentation)
  expect_equal(scores[1], 1)
  expect_equal(scores[2], 1)
})

test_that("reproductible results", {
  set.seed(10)
  load("test_data/complex_data.rda")
  true_segmentation <- rep(1:5, each = 20)
  result <- mlcc.bic(x, numb.clusters = 4:6, max.dim = 5, numb.cores = 1)
  set.seed(10)
  expect_equal(result$segmentation, mlcc.bic(x, numb.clusters = 4:6, max.dim = 5, numb.cores = 1)$segmentation)
})

test_that("casting data frame to matrix", {
  load("test_data/small_matrix.rda")
  X <- data.frame(X)
  expect_warning(mlcc.bic(X, numb.clusters = 1:1, max.dim = 2, numb.cores = 1), "X is not a matrix*")
})

test_that("missing values", {
  load("test_data/small_matrix.rda")
  X[10, 10] <- NaN
  expect_warning(mlcc.bic(X, numb.clusters = 1:1, max.dim = 2, numb.cores = 1), "Missing values are*")
})

test_that("non numeric values", {
  load("test_data/small_matrix.rda")
  X[10, 10] <- "nonnumeric"
  expect_error(mlcc.bic(X, numb.clusters = 1:1, max.dim = 2, numb.cores = 1), "*The following variables are not quantitative*")
})
