#' Variable Clustering with Multiple Latent Components Clustering algorithm
#'
#' Package varclust performs clustering of variables, according
#' to a probabilistic model, which assumes that each clusters lies in a low
#' dimensional subspace. Segmentation of variables, number of clusters and their
#' dimensions are selected based on the appropriate implementation of the 
#' Bayesian Information Criterion. 
#'
#' The best candidate models are identified by the specific 
#' implementation of K-means algorithm, in which cluster centers are represented
#' by some number of orthogonal factors and distance between a given variable
#' and a cluster center depends on residuals from a linear model fit. Based on
#' the Bayesian Information Criterion (BIC),  sums of squares of residuals are 
#' appropriately scaled, which allows to avoid an over-excessive attraction by
#' clusters with larger dimensions. To reduce the chance that the local minimum
#' of BIC is obtained instead of the global one, for every fixed number of 
#' clusters in a given range K-means algorithm is run large number of times,
#' with different random initializations of cluster centers.
#' 
#' The main function of package \pkg{varclust} is
#' \code{\link{mlcc.bic}} which allows clustering variables in a data 
#' with unknown number of clusters. Variable partition is computed
#' with k-means based algorithm. Number of clusters and their dimensions
#' are computed using BIC criterion.
#' If the number of clusters is known one might use function \code{\link{mlcc.reps}},
#' which takes number of clusters as a parameter. For \code{\link{mlcc.reps}} one might 
#' specify as well some initial segmentation for k-means algorithm. This can be useful if
#' user has some apriori knowledge about clustering.
#' 
#' We also provide function \code{\link{misclassification}} that computes misclassification
#' rate between two partitions. This performance measure is 
#' extensively used in image segmentation.
#'
#' @docType package
#' @name varclust
#' @details Version: 0.9.23
#' @importFrom RcppEigen fastLmPure
#' @importFrom doMC registerDoMC
#' @importFrom parallel detectCores
#' @import foreach
#' @author Piotr Sobczyk,
#'         Julie Josse,
#'         Malgorzata Bogdan
#' 
#' Maintainer: Piotr Sobczyk \email{Piotr.Sobczyk@@pwr.edu.pl}
#' 
#' @examples
#' \donttest{
#' sim.data <- data.simulation(n = 50, SNR = 1, K = 3, numb.vars = 50, max.dim = 3)
#' mlcc.bic(sim.data$X, numb.clusters = 1:5, numb.runs = 20)
#' mlcc.reps(sim.data$X, numb.clusters = 3, numb.runs = 20)}
NULL