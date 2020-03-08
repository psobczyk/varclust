#' Variable Clustering with Multiple Latent Components Clustering algorithm
#' 
#' Package varclust performs clustering of variables, according to a 
#' probabilistic model, which assumes that each cluster lies in a low 
#' dimensional subspace. Segmentation of variables, number of clusters and their
#' dimensions are selected based on the appropriate implementation of the 
#' Bayesian Information Criterion.
#' 
#' The best candidate models are identified by the specific implementation of 
#' K-means algorithm, in which cluster centers are represented by some number of
#' orthogonal factors(principal components of the variables within a cluster) 
#' and similarity between a given variable and a cluster center depends on 
#' residuals from a linear model fit. Based on the Bayesian Information 
#' Criterion (BIC), sums of squares of residuals are appropriately scaled, which
#' allows to avoid an over-excessive attraction by clusters with larger 
#' dimensions. To reduce the chance that the local minimum of modified BIC 
#' (mBIC) is obtained instead of the global one, for every fixed number of 
#' clusters in a given range K-means algorithm is run large number of times, 
#' with different random initializations of cluster centers.
#' 
#' The main function of package \pkg{varclust} is \code{\link{mlcc.bic}} which 
#' allows clustering variables in a data with unknown number of clusters. 
#' Variable partition is computed with k-means based algorithm. Number of 
#' clusters and their dimensions are estimated using mBIC and PESEL 
#' respectively. If the number of clusters is known one might use function 
#' \code{\link{mlcc.reps}}, which takes number of clusters as a parameter. For 
#' \code{\link{mlcc.reps}} one might specify as well some initial segmentation 
#' for k-means algorithm. This can be useful if user has some a priori knowledge
#' about clustering.
#' 
#' We provide also two functions to simulate datasets with described structure. 
#' The function \code{\link{data.simulation}} generates the data so that the 
#' subspaces are indepentend and \code{\link{data.simulation.factors}} generates
#' the data where some factores are shared between the subspaces.
#' 
#' We also provide function measures of quality of clustering. 
#' \code{\link{misclassification}} computes misclassification rate between two 
#' partitions. This performance measure is extensively used in image 
#' segmentation. The other measure is implemented as \code{\link{integration}} 
#' function.
#' 
#' @docType package
#' @name varclust
#' @details Version: 0.9.5
#' @importFrom RcppEigen fastLmPure
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom parallel detectCores
#' @importFrom pesel pesel
#' @importFrom graphics axis plot plot.default points
#' @importFrom stats cov dnorm pnorm prcomp rnorm runif var
#' @importFrom utils str
#' @importFrom Matrix bdiag
#' @import doRNG
#' @import foreach
#' @author Piotr Sobczyk, Stanislaw Wilczynski, Julie Josse, Malgorzata Bogdan
#'   
#'   Maintainer: Piotr Sobczyk \email{pj.sobczyk@@gmail.com}
#'   
#' @examples
#' \donttest{
#' sim.data <- data.simulation(n = 50, SNR = 1, K = 3, numb.vars = 50, max.dim = 3)
#' mlcc.bic(sim.data$X, numb.clusters = 1:5, numb.runs = 20, numb.cores = 1, verbose = TRUE)
#' mlcc.reps(sim.data$X, numb.clusters = 3, numb.runs = 20, numb.cores = 1)
#' }
NULL