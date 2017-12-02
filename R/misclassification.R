#' Computes misclassification rate
#' 
#' Missclasification is a commonly used performance measure in subspace clustering.
#' It allows to compare two partitions with the same number of clusters.
#'  
#' As getting exact value of misclassification requires checking all permutations 
#' and is therefore intrackable even for modest number of clusters, a heuristic approach is proposed.
#' It is assumed that there are K classes of maximum M elements. 
#' Additional requirement is that classes labels are from range [1, K].
#' 
#' @param group a vector, first partition
#' @param true_group a vector, second (reference) partition
#' @param M an integer, maximal number of elements in one class
#' @param K an integer, number of classes
#' @references {R. Vidal. Subspace clustering. Signal Processing Magazine, IEEE, 28(2):52-68,2011}
#' @export
#' @return misclassification rate
#' @examples
#' \donttest{
#' sim.data <- data.simulation(n = 100, SNR = 1, K = 5, numb.vars = 30, max.dim = 2)
#' mlcc.fit <- mlcc.reps(sim.data$X, numb.clusters = 5, numb.runs = 20, max.dim = 2)
#' misclassification(mlcc.fit$segmentation,sim.data$s, 30, 5)
#' }
#' 
#' #one can use this function not only for clusters
#' partition1 <- sample(10, 300, replace = TRUE)
#' partition2 <- sample(10, 300, replace = TRUE)
#' misclassification(partition1, partition1, max(table(partition1)), 10)
#' misclassification(partition1, partition2, max(table(partition2)), 10)
misclassification <-function(group, true_group, M, K){
  if (length(group) != length(true_group))
    stop("Partitions are of different lengths")
  forbidden <- NULL
  suma <- 0
  nG = max(group);
  for (i in M:1) { #differnet concordance levels
    for (j in 1:nG) { #subspaces numbers (found)
      if (sum(j==forbidden)==0) { #subspace not yet used
        for (k in 1:K) { # subspaces numbers (true)
          if (sum(j==group[true_group==k])==i) {
            suma <- suma + i
            forbidden <- c(forbidden, j)
            break
          }
        }
      }
    }
  }
  mis <- 1-suma/length(true_group)
  return(mis)
}


#' Computes integration of each cluster in given partition
#' 
#' 
#' @param group a vector, first partition
#' @param true_group a vector, second (reference) partition
#' @references {Praca Sołtysa}
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