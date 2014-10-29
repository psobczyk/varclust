#' Variable Clustering with Multiple Latent Components Clustering algorithm
#'
#' Variable Clustering with Multiple Latent Components Clustering is based on k-means algorithm. 
#' In each step cluster centers are few PCA components, computed for variables in that cluster. 
#' Than distance is defined by R^2 (obtained by performing least-squares).
#'
#' The main function of package \pkg{VarClust} is
#' \code{\link{mlcc.bic}} which allows you to cluster data 
#' with unknown number of clusters. Variable partition is computed
#' with k-means based algorithm. Number of clusters and their dimensions
#' are computed using BIC criterion.
#' If user knows the number of clusters one might use function \code{\link{mlcc.reps}},
#' which takes number of clusters as a parameter. For \code{\link{mlcc.reps}} one might 
#' specify as well some initial segmentation for k-means algorithm. This can be useful if
#' user has some apriori knowledge about clustering.
#' 
#' We also provide function \code{\link{misclassification}} that computes misclassification
#' rate between two partitions in which might have different labelings. This performance measure is 
#' extensively used in application to computer vision.
#'
#' @docType package
#' @name VarClust
#' @details Version: 0.9.8
#' @importFrom RcppEigen fastLmPure
#' @importFrom doMC registerDoMC
#' @importFrom parallel detectCores
#' @import foreach
#' @author Piotr Sobczyk,
#'         Julie Josse
#' 
#' Maintainer: Piotr Sobczyk \email{Piotr.Sobczyk@@pwr.edu.pl}
#' @references 
#' Malgorzata Bogdan, Julie Josse, Piotr Sobczyk, \emph{Clustering around latent variables - a technical report}, 2014, 
#' \url{www.im.pwr.edu.pl/~sobczyk/research.html} 
#' 
#' @examples
#' \donttest{
#' data <- data.simulation(n=100, SNR=1, K=5, numbVars=30, max.dim=2)
#' mlcc.bic(data$X, numb.clusters=1:5, numb.runs=20)
#' mlcc.reps(data$X, numb.clusters=5, numb.runs=20, max.dim=3)}
NULL