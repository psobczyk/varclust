##' @title Computing joint sigma for all clusters
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

#' Compute missclasification rate for subspace clustering
#' 
#' As getting exact value requires checking all permutations a heuristic approach is proposed
#' It is assumed that there are n cluster each of N elements and that they are sorted.
#' 
#' @param group proposed clustering
#' @param N number of elements in each cluster
#' @param n number of clusters
#' @export
#' @return mis misclassification rate
missclassify.heuristic <-function(group, N, n){
  forbidden = NULL;
  suma = 0;
  nG = max(group);
  for (i in N:1){ #differnet concordance levels
    for(j in 1:nG){ #subspaces numbers (found)
      if (sum(j==forbidden)==0){ #subspace not yet used
        for (k in 1:n){ # subspaces numbers (true)
          if (sum(j==group[((k-1)*N+1):(k*N)])==i){
            suma = suma + i
            forbidden = c(forbidden, j)
            break;
          }
        }
      }
    }
  }
  mis = 1-suma/(N*n)
  return(mis)
}

#' Computes missclasification rate.
#' 
#' As getting exact value of missclasification requires checking all permutations 
#' and is therefore intrackable even for modest number of classes, a heuristic approach is proposed.
#' It is assumed that there are n classes of maximum N elements. 
#' Additional requirement is that classes identifiers are from range [1, n]
#' 
#' @param group a vector, first partiton
#' @param true_group a vector, second (reference) partition
#' @param N an integer, maximal number of elements in one class
#' @param n an integer, number of classes
#' @export
#' @return misclassification rate
#' @examples
#' partition1 <- sample(10, 300, replace=T)
#' partition2 <- sample(10, 300, replace=T)
#' missclassify.heuristic2(partition1, partition1, max(table(partition1)), 10)
#' missclassify.heuristic2(partition1, partition2, max(table(partition2)), 10)
missclassify.heuristic2 <-function(group, true_group, N, n){
  forbidden = NULL;
  suma = 0;
  nG = max(group);
  for (i in N:1){ #differnet concordance levels
    for(j in 1:nG){ #subspaces numbers (found)
      if (sum(j==forbidden)==0){ #subspace not yet used
        for (k in 1:n){ # subspaces numbers (true)
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
