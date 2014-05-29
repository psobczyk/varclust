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