#' Computes integration of each cluster in given partition
#' 
#' @param group a vector, first partition
#' @param true_group a vector, second (reference) partition
#' @references {M. Sołtys. Metody analizy skupień. Master’s thesis, Wrocław University of Technology, 2010}
#' @export
#' @return integration
#' 
integration <- function(group, true_group){
  n <- length(group)
  K <- max(unique(true_group))
  if (max(unique(group)) > K){
    stop("Number of clusters in given partition has to be less or equal to the number of clusters in true partition")
  }
  if (n != length(true_group))
    stop("Partitions are of different lengths")
  integrationMatrix <- matrix(0,nrow = K, ncol = K)
  for (i in 1:n){
    integrationMatrix[group[i],true_group[i]] = integrationMatrix[group[i],true_group[i]] + 1
  }
  clusters <- apply(integrationMatrix, 2, max)
  cluster_indices <- apply(integrationMatrix, 2, which.max)
  sizes_true <- apply(integrationMatrix, 2, sum)
  sizes_group <- apply(integrationMatrix, 1, sum) 
  int <- clusters/sizes_true
  acont <- clusters/sizes_group[cluster_indices]
  return(c(mean(int),mean(acont)))
}
