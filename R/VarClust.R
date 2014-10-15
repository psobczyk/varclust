#' Variable Clustering
#'
#' Variable Clustering based on k-means algorithm. In each step cluster center is a number of PCA components
#' of points in that cluster. Than distance is defined as R^2 in least-squares.
#'
#' The main function of package \pkg{VarClust} is
#' \code{\link{MPCV.BIC}} which allows you to cluster data 
#' with unknown number of clusters. If you have know
#' the number of clusters use function \code{\link{MPCV.reps}}
#'
#' @docType package
#' @name VarClust
#' @details Version: 0.9
#' @importFrom RcppEigen, mclust, ClustOfVar, foreach, doMC
#' @author Piotr Sobczyk
#' 
#' Maintainer: Piotr Sobczyk \email{Piotr.Sobczyk@@pwr.edu.pl}
#' @references Ulrike von Luxburg 2007 \emph{A Tutorial on Spectral Clustering} 
#' 
#' Marie Chavent, Benoit Liquet, Vanessa Kuentz-Simonet, Jerome Saracco 2013 \emph{ClustOfVar: An R Package for the Clustering of Variables} 
NULL