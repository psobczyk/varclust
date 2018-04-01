context("Testing mlcc.bic")

library(varclust)

test_that("mlcc.bic for small matrix",{
  load("test_data/small_matrix.rda")
  true_segmentation <- rep(1:2, each=50)
  result <- mlcc.bic(X, numb.clusters = 1:5, max.dim = 2, numb.cores = 1)
  expect_equal(2, result$nClusters)
  scores <- varclust::integration(result$segmentation, true_segmentation)
  expect_gte(scores[1],0.9)
  expect_gte(scores[2],0.9)
})
