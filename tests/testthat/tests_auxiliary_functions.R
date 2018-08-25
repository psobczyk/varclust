context("Testing auxiliary functions")

library(varclust)

test_that("missclassification works ok", {
  part1 <- c(rep(1, 10), rep(2,10))
  part2 <- c(rep(1,15), rep(2,5))
  expect_equal(misclassification(part1, part2, 15, 2), 0.25)
  expect_equal(misclassification(part1, part2, 15, 2), misclassification(part2, part1, 15, 2))
  expect_equal(misclassification(part1, part2, 25, 2), misclassification(part2, part1, 15, 2))
  expect_false(misclassification(part1, part2, 8, 2) == misclassification(part2, part1, 15, 2))
  expect_error(misclassification(part1, c(3, part2), 15, 2))
})
test_that("integration",{
  part1 <- c(rep(1, 10), rep(2,10))
  part2 <- c(rep(1,15), rep(2,5))
  expect_equal(integration(part2,part1)[1], 0.75 )
  expect_equal(integration(part2,part1)[2], 0.5 )
})

test_that("choose cluster BIC", {
  set.seed(1)
  sim.data <- data.simulation(n = 100, SNR = 1, K = 2, numb.vars = 4, max.dim = 2)
  pca1 <- summary(prcomp(x=sim.data$signals[,1:4]))$x[,1:2]
  pca2 <- summary(prcomp(x=sim.data$signals[,5:8]))$x[,1:2]
  expect_equal(choose.cluster.BIC(sim.data$signals[,1], list(pca1, pca2), 2), 1)
  expect_equal(choose.cluster.BIC(sim.data$signals[,5], list(pca1, pca2), 2), 2)
})

test_that("choose cluster BIC - when n is large we can choose the right model", {
  set.seed(1)
  sim.data <- data.simulation(n = 100, SNR = 10, K = 2, numb.vars = 10, max.dim = 2)
  pca1 <- summary(prcomp(x=sim.data$X[,1:10]))$x[,1:2]
  pca2 <- summary(prcomp(x=sim.data$X[,11:20]))$x[,1:2]
  expect_equal(choose.cluster.BIC(sim.data$X[,1], list(pca1, pca2), 2), 1)
  expect_equal(choose.cluster.BIC(sim.data$X[,11], list(pca1, pca2), 2), 2)
})

 test_that("choose cluster BIC - warning for suspiciously low noise", {
   set.seed(1)
   sim.data <- data.simulation(n = 100, SNR = 1, K = 2, numb.vars = 4, max.dim = 2)
   pca1 <- summary(prcomp(x=sim.data$signals[,1:4]))$x[,1:2]
   pca2 <- summary(prcomp(x=sim.data$signals[,1:4]))$x[,1:3]
   expect_warning(choose.cluster.BIC(sim.data$signals[,1], list(pca1, pca2), 2, TRUE))
 })

test_that("sim.data - factors are orthogonal", {
  set.seed(1)
  sim.data <- data.simulation(n = 50, SNR = 1, K = 2, numb.vars = 20, max.dim = 2)
  expect_gt(cor.test(sim.data$factors[,1], sim.data$factors[,2])$p.value, 0.05)
  expect_gt(cor.test(sim.data$factors[,3], sim.data$factors[,4])$p.value, 0.05)
})

test_that("warning for -inf", {
  load("test_data/small_matrix.rda")
  segmentation <- c(rep(1, each=98),2,2)
  expect_warning(cluster.pca.BIC(X, segmentation, c(5,3), 2, 7), "The dimensionality *")
})
