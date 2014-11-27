#' Multiple Latent Components Clustering - kmeans algorithm
#'
#' Performs k-means based subspace clustering. Center of each cluster is some number 
#' of principal components. Similarity measure is R^2 coefficient.
#'
#' @param X a matrix with only continuous variables
#' @param number.clusters an integer, number of clusters to be used
#' @param stop.criterion an integer indicating how many changes in partitions triggers stopping the algorithm
#' @param max.iter an integer, maximum number of iterations of k-means
#' @param max.subspace.dim an integer, maximum dimension of subspaces
#' @param initial.segmentation a vector, initial segmentation of variables to clusters
#' @param estimate.dimensions a boolean, if TRUE (value set by default) subspaces dimensions are estimated
#' @export
#' @return A list consisting of:
#' \item{segmentation}{a vector containing the partition of the variables}
#' \item{pcas}{a list of matrices, basis vectors for each cluster (subspace)}
mlcc.kmeans <- function(X, number.clusters=2, stop.criterion=1, max.iter=40, max.subspace.dim=4, initial.segmentation=NULL, estimate.dimensions=FALSE){
  numbVars = dim(X)[2]
  rowNumb = dim(X)[1]
  pcas <- list(NULL)
  
  if(is.null(initial.segmentation)){
    los = sample(1:numbVars,number.clusters)
    pcas = lapply(1:number.clusters, function(i) matrix(X[,los[i]], nrow=rowNumb))
    segmentation <- sapply(1:numbVars,  function(j) choose.cluster(X[,j],pcas, number.clusters))
  }
  else
    segmentation=initial.segmentation
  
  new.segmentation <- segmentation
  for (iter in 1:max.iter){
    pcas <- lapply(1:number.clusters, function(i){
      if(dim(X[,segmentation==i, drop=F])[2]>max.subspace.dim){
        a <- summary(prcomp(x=X[,segmentation==i]))
        if (estimate.dimensions) {
          cut <- which.max(sapply(1:max.subspace.dim, 
                                  function(dim) cluster.BIC(X[,segmentation==i], rep(1,sum(segmentation==i)), max.dim = dim, numb.clusters = 1)))
        }
        else {
          cut <- max.subspace.dim
        }
        return(matrix(a$x[,1:cut], nrow=rowNumb))
      }
      else{
          return(matrix(rnorm(rowNumb*max.subspace.dim), nrow=rowNumb))
      }
    })
    if (estimate.dimensions)
      new.segmentation <- sapply(1:numbVars, function(j) choose.cluster.BIC(X[,j], pcas, number.clusters))
    else {
      new.segmentation <- sapply(1:numbVars, function(j) choose.cluster(X[,j], pcas, number.clusters))
    }
    if(sum(new.segmentation!=segmentation)<stop.criterion) break
    segmentation = new.segmentation
  }
  return(list(segmentation=segmentation, 
              pcas=pcas))
}