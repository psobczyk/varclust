#' Variable Clustering with Multiple Latent Components Clustering algorithm
#'
#' Variable Clustering with Multiple Latent Components Clustering is based on k-means algorithm. 
#' In each step cluster centers are few PCA components, computed for variables in that cluster. 
#' The distance is defined by R^2 (obtained by performing least-squares).
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
#' @details Version: 0.9.21
#' @importFrom RcppEigen fastLmPure
#' @importFrom doMC registerDoMC
#' @importFrom parallel detectCores
#' @import foreach
#' @author Piotr Sobczyk,
#'         Julie Josse
#' 
#' Maintainer: Piotr Sobczyk \email{Piotr.Sobczyk@@pwr.edu.pl}
#' @references 
#' Piotr Sobczyk, Malgorzata Bogdan, Julie Josse, \emph{Clustering around latent variables - a technical report}, 2014, 
#' \url{www.im.pwr.edu.pl/~sobczyk/research.html} 
#' 
#' @examples
#' \donttest{
#' sim.data <- data.simulation(n=100, SNR=1, K=5, numb.vars=30, max.dim=2)
#' mlcc.bic(sim.data$X, numb.clusters=1:5, numb.runs=20)
#' mlcc.reps(sim.data$X, numb.clusters=5, numb.runs=20)}
NULL