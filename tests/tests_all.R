library(testthat)
test_that('pierwszy test', {
  expect_equal(10, 10)
})

test_that('subspace are in fact low rank', {
  w <- dataSIMULATION(n = 100, SNR=1, K = 2, numbPoints = 50, max.dim = 2)
  expect_that(which.max(summary(prcomp(w$signals[,1:50]))$importance[3,]==1), is_equivalent_to(2))
  expect_that(which.max(summary(prcomp(w$signals[,51:100]))$importance[3,]==1), is_equivalent_to(2))
})

test_that('noisy subspace is not low rank', {
  w <- dataSIMULATION(n = 100, SNR=1, K = 2, numbPoints = 50, max.dim = 2)
  expect_true(which.max(summary(prcomp(w$X[,1:50]))$importance[3,]==1)>2)
  expect_true(which.max(summary(prcomp(w$X[,51:100]))$importance[3,]==1)>2)
})

test_that('proper choice of closest cluster', {
  n   <- 100
  pca <- NULL
  w <- dataSIMULATION(n = 100, SNR=1, K = 2, numbPoints = 50, max.dim = 2)
  pca[[1]] <- matrix(summary(prcomp(w$signal[,1:50]))$x[,1:2], nrow=n)
  pca[[2]] <- matrix(summary(prcomp(w$signal[,51:100]))$x[,1:2], nrow=n)
  expect_true(choose.cluster(w$signal[,1], pca, 2)==1) 
  expect_false(choose.cluster(w$signal[,23], pca, 2)==2)
  expect_true(choose.cluster(w$signal[,51], pca, 2)==2)
  expect_false(choose.cluster(w$signal[,51], pca, 2)==1)
})
