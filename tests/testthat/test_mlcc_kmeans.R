context("Testing mlcc.kmeans")

library(varclust)

test_that("principal components returned by mlcc.kmeans are in fact principal components of the clusters", {
  X <- varclust::data.simulation.factors(n=20, K = 2, numb.vars = 50, numb.factors = 5)$X
  result <- mlcc.kmeans(X, number.clusters = 2, max.iter = 1)
  segmentation <- result[[1]]
  pcas <- result[[2]]
  X1 <- X[,segmentation==1, drop=F]
  X2 <- X[,segmentation==2, drop=F]
  dims <- sapply(pcas, ncol)
  real_pcas1 <- matrix(summary(prcomp(x=X1))$x[,1:dims[1]], nrow = 20)
  real_pcas2 <- matrix(summary(prcomp(x=X2))$x[,1:dims[2]], nrow = 20)
  expect_equal(real_pcas1, pcas[[1]])
  expect_equal(real_pcas2, pcas[[2]])
  
})
