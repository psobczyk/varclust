context("Testing mlcc.reps")

library(varclust)

test_that("mlcc.reps on small matrix", {
  set.seed(1)
  load("test_data/small_matrix.rda")
  result <- varclust::mlcc.reps(X, max.dim = 2, numb.cores = 1)
  true_segmentation <- rep(1:2, each=50)
  scores <- varclust::integration(result$segmentation, true_segmentation)
  expect_gte(scores[1],0.9)
  expect_gte(scores[2],0.9)
})
