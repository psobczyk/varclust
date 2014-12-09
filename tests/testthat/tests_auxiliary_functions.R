context("Testing auxiliary functions")

library(varclust)

test_that("missclassification works ok", {
  part1 <- c(rep(1, 10), rep(2,10))
  part2 <- c(rep(1,15), rep(2,5))
  expect_equal(varclust::misclassification(part1, part2, 15, 2), 0.25)
  expect_equal(varclust::misclassification(part1, part2, 15, 2), varclust::misclassification(part2, part1, 15, 2))
  expect_equal(varclust::misclassification(part1, part2, 25, 2), varclust::misclassification(part2, part1, 15, 2))
  expect_false(varclust::misclassification(part1, part2, 8, 2) == varclust::misclassification(part2, part1, 15, 2))
  expect_error(varclust::misclassification(part1, c(3, part2), 15, 2))
})

test_that("choose cluster", {
  sim.data <- varclust::data.simulation(n = 100, SNR = 1, K = 2, numb.vars = 4, max.dim = 2)
  pca1 <- summary(prcomp(x=sim.data$signals[,1:4]))$x[,1:2]
  pca2 <- summary(prcomp(x=sim.data$signals[,5:8]))$x[,1:2]
  expect_equal(varclust:::choose.cluster(sim.data$signals[,1], list(pca1, pca2), 2), 1)
  expect_equal(varclust:::choose.cluster(sim.data$signals[,5], list(pca1, pca2), 2), 2)
})

test_that("choose cluster BIC - when n is large we can choose the right model", {
  sim.data <- varclust::data.simulation(n = 100, SNR = 10, K = 2, numb.vars = 10, max.dim = 2)
  pca1 <- summary(prcomp(x=sim.data$X[,1:10]))$x[,1:2]
  pca2 <- summary(prcomp(x=sim.data$X[,11:20]))$x[,1:2]
  expect_equal(varclust:::choose.cluster.BIC(sim.data$X[,1], list(pca1, pca2), 2), 1)
  expect_equal(varclust:::choose.cluster(sim.data$X[,11], list(pca1, pca2), 2), 2)
})

test_that("choose cluster BIC - warning for suspiciously low noise", {
  sim.data <- varclust::data.simulation(n = 100, SNR = 1, K = 2, numb.vars = 4, max.dim = 2)
  pca1 <- summary(prcomp(x=sim.data$signals[,1:4]))$x[,1:2]
  pca2 <- summary(prcomp(x=sim.data$signals[,1:4]))$x[,1:3]
  expect_warning(varclust:::choose.cluster.BIC(sim.data$signals[,1], list(pca1, pca2), 2))
})


test_that("adjusted cluster BIC - is computed correctly. Works with mlcc.reps", {
  set.seed(1)
  sim.data <- varclust::data.simulation(n = 50, SNR = 1, K = 2, numb.vars = 20, max.dim = 2)
  expect_equal(varclust:::adjusted.cluster.BIC(scale(sim.data$X), sim.data$s, c(2,2), 2), 
               -2457.157, tolerance = 1e-3, scale = 1)
  #please note that for the mlcc.fit value is -2445.326
  #the difference is because of sigma estimation. In mlcc.fit we get smaller sigma because
  #we estimate it with max.dim, hence BIC is higher.
})

test_that("get sigma", {
  set.seed(1)
  sim.data <- varclust::data.simulation(n = 50, SNR = 1, K = 2, numb.vars = 20, max.dim = 2)
  sigma2 <- varclust:::getSigma(X = scale(sim.data$X), segmentation = sim.data$s, max.dim = 2, 
                     n = 50, p = 2*20, numb.clusters = 2)
  sigma3 <- varclust:::getSigma(X = scale(sim.data$X), segmentation = sim.data$s, max.dim = 3, 
                     n = 50, p = 2*20, numb.clusters = 2)
  expect_more_than(sigma2, sigma3)
  sigma2_sig <- varclust:::getSigma(X = scale(sim.data$signals), segmentation = sim.data$s, max.dim = 2, 
                               n = 50, p = 2*20, numb.clusters = 2)
  expect_equal(sigma2_sig, 0)
  # the values below were computed using version 0.9.19
  expect_equal(sigma2, 0.7088557, tolerance = 1e-7, scale = 1)
  expect_equal(sigma3, 0.6795928, tolerance = 1e-7, scale = 1)
})

test_that("sim.data - factors are orthogonal", {
  set.seed(1)
  sim.data <- varclust::data.simulation(n = 50, SNR = 1, K = 2, numb.vars = 20, max.dim = 2)
  expect_more_than(cor.test(sim.data$factors[,1], sim.data$factors[,2])$p.value, 0.05)
  expect_more_than(cor.test(sim.data$factors[,3], sim.data$factors[,4])$p.value, 0.05)
})