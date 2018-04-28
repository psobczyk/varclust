#' Computes integration and acontamination of the clustering
#' 
#' Integartion and acontamination are measures of the quality of a clustering 
#' with a reference to a true partition. Let \eqn{X = (x_1, \ldots x_p)} be the 
#' data set, \eqn{A} be a partition into clusters \eqn{A_1, \ldots A_n} (true 
#' partition) and \eqn{B} be a partition into clusters \eqn{B_1, \ldots, B_m}. 
#' Then for cluster \eqn{A_j} integration is eqaul to: \deqn{Int(A_j) = 
#' \frac{max_{k = 1, \ldots, m} \# \{  i \in \{ 1, \ldots p \}: x_i \in A_j 
#' \wedge x_i \in B_k \}  }{\# A_j}} The \eqn{B_k} for which the value is 
#' maximized is called the integrating cluster of \eqn{A_j}. Then the 
#' integration for the whole clustering equals is \eqn{Int(A,B) = \frac{1}{n} 
#' \sum_{j=1}^n Int(A_j)} .The acontamination is defined by: \deqn{Acont(A_j) = 
#' \frac{ \# \{  i \in \{ 1, \ldots p \}: x_i \in A_j \wedge x_i \in B_k \} }{\#
#' B_k}} where \eqn{B_k} is the integrating cluster for \eqn{A_j}. Then the 
#' acontamination for the whole dataset is \eqn{Acont(A,B) = \frac{1}{n} 
#' \sum_{j=1}^n Acont(A_j)}
#' 
#' @param group A vector, first partition.
#' @param true_group A vector, second (reference) partition.
#' @references {M. Sołtys. Metody analizy skupień. Master’s thesis, Wrocław 
#'   University of Technology, 2010}
#' @export
#' @return An array containing values of integration and acontamination.
#' @examples
#' \donttest{
#' sim.data <- data.simulation(n = 20, SNR = 1, K = 2, numb.vars = 50, max.dim = 2)
#' true_segmentation <- rep(1:2, each=50)
#' mlcc.fit <- mlcc.reps(sim.data$X, numb.clusters = 2, max.dim = 2, numb.cores=1)
#' integration(mlcc.fit$segmentation, true_segmentation)}
#' 
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
