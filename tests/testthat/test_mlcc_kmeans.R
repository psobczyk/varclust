context("Testing mlcc.kmeans")

library(varclust)

test_that("when data is random we select dimension equal to 0", {
  X <- with(set.seed(23), matrix(rnorm(1000), ncol=20))
  mlcc.res <- mlcc.kmeans(X, number.clusters = 2,
                          max.iter = 20, max.subspace.dim = 3)
  expect_equal(unique(mlcc.res$segmentation), 1)
  expect_equal(ncol(mlcc.res$pcas[[mlcc.res$segmentation[1]]]), 0)
})


test_that("principal components returned by mlcc.kmeans are in fact principal components of the clusters", {
  set.seed(10)
  X <- data.simulation.factors(n = 20, K = 2, numb.vars = 50, numb.factors = 5)$X
  result <- mlcc.kmeans(X, number.clusters = 2, max.iter = 1)
  segmentation <- result[[1]]
  pcas <- result[[2]]
  X1 <- X[, segmentation == 1, drop = F]
  X2 <- X[, segmentation == 2, drop = F]
  dims <- sapply(pcas, ncol)
  real_pcas1 <- matrix(summary(prcomp(x = X1))$x[, 1:dims[1]], nrow = 20)
  real_pcas2 <- matrix(summary(prcomp(x = X2))$x[, 1:dims[2]], nrow = 20)
  expect_equal(real_pcas1, pcas[[1]])
  expect_equal(real_pcas2, pcas[[2]])
})


test_that("incorrect length on initial segmentation", {
  load("test_data/small_matrix.rda")
  segmentation <- rep(1:2, each = 40)
  expect_error(mlcc.kmeans(X, number.clusters = 2, initial.segmentation = segmentation), "The lenght of initial*")
})


test_that("incorrect initial segmentation", {
  load("test_data/small_matrix.rda")
  segmentation <- rep(1:4, each = 25)
  expect_error(mlcc.kmeans(X, number.clusters = 2, initial.segmentation = segmentation), "Too many*")
})


test_that("perfect segmentation", {
  load("test_data/small_matrix.rda")
  true_segmentation <- rep(1:2, each = 50)
  segmentation <- mlcc.kmeans(X, initial.segmentation = true_segmentation)$segmentation
  expect_equal(segmentation, true_segmentation)
})
