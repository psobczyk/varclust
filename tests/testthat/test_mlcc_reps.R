context("Testing mlcc.reps")

library(varclust)

test_that("mlcc.reps on small matrix", {
  set.seed(1)
  load("test_data/small_matrix.rda")
  result <- mlcc.reps(X, max.dim = 2, numb.cores = 1)
  true_segmentation <- rep(1:2, each=50)
  scores <- integration(result$segmentation, true_segmentation)
  expect_gte(scores[1],0.9)
  expect_gte(scores[2],0.9)
})

test_that("perfect segmentation", {
  set.seed(1)
  load("test_data/small_matrix.rda")
  true_segmentation <- rep(1:2, each=50)
  bad_segmentation <- c(rep(1:2, each=25), rep(1:2, each=25))
  segmentation <- mlcc.reps(X, initial.segmentations = list(true_segmentation, bad_segmentation), numb.cores = 1)$segmentation
  expect_equal(segmentation, true_segmentation)
})


test_that("casting data frame to matrix", {
  load("test_data/small_matrix.rda")
  X <- data.frame(X)
  expect_warning(mlcc.reps(X, numb.clusters = 1, max.dim = 2, numb.cores = 1), "X is not a matrix*")
})

test_that("missing values", {
  load("test_data/small_matrix.rda")
  X[10,10] <- NaN
  expect_warning(mlcc.reps(X, numb.clusters = 1, max.dim = 2, numb.cores = 1), "Missing values are*")
})

test_that("non numeric values", {
  load("test_data/small_matrix.rda")
  X[10,10] <- "nonnumeric"
  expect_error(mlcc.reps(X, numb.clusters = 1:1, max.dim = 2, numb.cores = 1), "*The following variables are not quantitative*")
})
