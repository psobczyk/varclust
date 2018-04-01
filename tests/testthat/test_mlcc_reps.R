context("Testing mlcc.reps")

library(varclust)

test_that("mlcc.reps on small matrix", {
  load("test_data/small_matrix.rda")
  result <- varclust::mlcc.reps(X, max.dim = 2)
  true_segmentation <- rep(1:2, each=50)
  scores <- varclust::integration(result$segmentation, true_segmentation)
  print(scores)
  expect_gte(scores[1],0.9)
  expect_gte(scores[2],0.9)
})
