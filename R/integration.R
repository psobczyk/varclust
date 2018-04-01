#' Computes integration and acontamination of each cluster in given partition
#' 
#' @param group a vector, first partition
#' @param true_group a vector, second (reference) partition
#' @references {M. Sołtys. Metody analizy skupień. Master’s thesis, Wrocław University of Technology, 2010}
#' @export
#' @return (integration, acontamination)
#' @examples
#' \donttest{
#' sim.data <- data.simulation(n = 20, SNR = 1, K = 2, numb.vars = 50, max.dim = 2)
#' true_segmentation <- rep(1:2, each=50)
#' mlcc.fit <- mlcc.reps(sim.data$X, numb.clusters = 2, max.dim = 2)
#' integration(true_segmentation, mlcc.fit$segmentation)
#' }
integration <- function(group, true_group){
  n <- length(group)
  K1 <- max(unique(group))
  K2 <- max(unique(true_group))
  if (n != length(true_group))
    stop("Partitions are of different lengths")
  integrationMatrix <- matrix(0,nrow = K1, ncol = K2)
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
