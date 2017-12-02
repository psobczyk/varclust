#' Computes integration of each cluster in given partition
#' 
#' 
#' @param group a vector, first partition
#' @param true_group a vector, second (reference) partition
#' @references {TO DO: add reference to Michal Soltys doctoral dissertation}
#' @export
#' @return integration
#' 
integration <- function(group, true_group){
  n <- length(group)
  K <- max(unique(true_group))
  if (n != length(true_group)) # || !all(sort(unique(group)) == sort(unique(true_group)))) - zdarzają się puste grupy
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
