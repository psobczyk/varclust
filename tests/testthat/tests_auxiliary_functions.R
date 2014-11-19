context("Testing auxiliary functions")

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
  expect_equal(varclust::choose.cluster(sim.data$signals[,1], list(pca1, pca2), 2), 1)
  expect_equal(varclust::choose.cluster(sim.data$signals[,5], list(pca1, pca2), 2), 2)
})

test_that("choose cluster BIC - when n is large we can choose the right model", {
  sim.data <- varclust::data.simulation(n = 100, SNR = 10, K = 2, numb.vars = 10, max.dim = 2)
  pca1 <- summary(prcomp(x=sim.data$X[,1:10]))$x[,1:2]
  pca2 <- summary(prcomp(x=sim.data$X[,11:20]))$x[,1:2]
  expect_equal(varclust::choose.cluster.BIC(sim.data$X[,1], list(pca1, pca2), 2), 1)
  expect_equal(varclust::choose.cluster(sim.data$X[,11], list(pca1, pca2), 2), 2)
})

test_that("choose cluster BIC - warning for suspiciously low noise", {
  sim.data <- varclust::data.simulation(n = 100, SNR = 1, K = 2, numb.vars = 4, max.dim = 2)
  pca1 <- summary(prcomp(x=sim.data$signals[,1:4]))$x[,1:2]
  pca2 <- summary(prcomp(x=sim.data$signals[,1:4]))$x[,1:3]
  expect_warning(varclust::choose.cluster.BIC(sim.data$signals[,1], list(pca1, pca2), 2))
})


test_that("cluster BIC - is computed correctly. Works with mlcc.reps", {
  set.seed(1)
  sim.data <- varclust::data.simulation(n = 50, SNR = 1, K = 2, numb.vars = 20, max.dim = 2)
  expect_true(abs(varclust::adjusted.cluster.BIC(scale(sim.data$X), sim.data$s, c(2,2), 2) == -2457.157)<1e-3)
  #please note that for the mlcc.fit value is -2445.326
  #the difference is because of sigma estimation. In mlcc.fit we get smaller sigma because
  #we estimate it with max.dim, hence BIC is higher.
})

