##' Computing joint sigma for all clusters
##'
##' Computes unbiased noise estimator under assumption that 
##' all subspaces are of the same dimension
##'
##' @param X data
##' @param segmentation variable segmentation
##' @param max.dim maximal subspace dimension
##' @param n number of rows
##' @param p number of variables
##' @param numb.clusters total number of clusters
##' @return unbiased noise estimator
##' @author Piotr Sobczyk
##' @keywords internal
getSigma <- function(X, segmentation, max.dim, n, p, numb.clusters){
  RES.sigma=0
  for(k in 1:numb.clusters){
    Xk = X[,segmentation==k, drop=F]
    if(dim(Xk)[2]>max.dim){ #length because it might be onedimensional
      svdSIGNAL= svd(Xk)  
      SIGNAL = matrix(svdSIGNAL$u[, 1:max.dim], ncol=max.dim) %*% 
        diag(svdSIGNAL$d[1:max.dim], nrow=max.dim) %*% 
        t(matrix(svdSIGNAL$v[, 1:max.dim], ncol=max.dim))
      RES.sigma = RES.sigma + sum((Xk - SIGNAL)^2)
    }
  }
  degrees.freedom <- n*p-p-n*max.dim-p*max.dim+max.dim^2+max.dim
  sigma <- sqrt(RES.sigma/degrees.freedom)
  sigma
}


#' Choose subspace closest to the given variable
#' 
#' The most similar subspace is choosen based on R^2
#' 
#' @param variable variable to be assigned
#' @param pcas orthogonal basis for different subspaces
#' @param numberClusters number of subspaces (clusters)
#' @return index number of subspace closest to variable
#' @keywords internal
choose.cluster <- function(variable, pcas, numberClusters){
  v1 = var(variable)
  which.max( vapply(1:numberClusters, function(i){
    v2 <- var(fastLmPure(pcas[[i]], variable, method = 0L)$residuals);
    p <- ncol(pcas[[i]]); 
    n <- length(variable);
    (v1-v2)/v1
  }, 0.9) )
}

#' Computes misclassification rate.
#' 
#' Missclasification is commonly used performance measure in subspace clustering.
#' It allows to compare two partitions with the same number of clusters.
#'  
#' As getting exact value of misclassification requires checking all permutations 
#' and is therefore intrackable even for modest number of clusters, a heuristic approach is proposed.
#' It is assumed that there are K classes of maximum M elements. 
#' Additional requirement is that classes labels are from range [1, K].
#' 
#' @param group a vector, first partiton
#' @param true_group a vector, second (reference) partition
#' @param M an integer, maximal number of elements in one class
#' @param K an integer, number of classes
#' @references {R Vidal. Subspace clustering. Signal Processing Magazine, IEEE, 28(2):52-68,2011.}
#' @export
#' @return misclassification rate
#' @examples
#' \donttest{
#' data <- dataSIMULATION(n=100, SNR=1, K=5, numbVars=30, max.dim=2)
#' mlcc.fit <- MPCV.reps(data$X, numb.clusters=5, numb.runs=20, max.dim=2)
#' misclassification(mlcc.fit$segmentation,data$s, 30, 5)
#' }
#' 
#' #one can use this function not only for clusters
#' partition1 <- sample(10, 300, replace=TRUE)
#' partition2 <- sample(10, 300, replace=TRUE)
#' misclassification(partition1, partition1, max(table(partition1)), 10)
#' misclassification(partition1, partition2, max(table(partition2)), 10)
misclassification <-function(group, true_group, M, K){
  forbidden = NULL;
  suma = 0;
  nG = max(group);
  for (i in M:1){ #differnet concordance levels
    for(j in 1:nG){ #subspaces numbers (found)
      if (sum(j==forbidden)==0){ #subspace not yet used
        for (k in 1:K){ # subspaces numbers (true)
          if (sum(j==group[true_group==k])==i){
            suma = suma + i
            forbidden = c(forbidden, j)
            break;
          }
        }
      }
    }
  }
  mis = 1-suma/length(true_group)
  return(mis)
}


#' Plot mlcc.fit class object
#' 
#' @param x mlcc.fit class object
#' @export
#' @keywords internal
plot.mlcc.fit <- function(x,...){
  clusterNumbs <- lapply(x$all.fit, function(y) y$nClusters)
  BICVals <- lapply(x$all.fit, function(y) y$BIC)
  plot.default(clusterNumbs, BICVals, type="b", xaxt="n", ylab="BIC", xlab="Number of clusters")
  axis(side = 1, labels = clusterNumbs, at=clusterNumbs)
}
