#' BIC for subspace clustering
#' 
#' Computes the value of BIC criterion for given data set and partition.
#' In each cluster we assume that variables are spanned by few factors.
#' Considering maximum likelihood we get that those factors are in fact
#' principal components. Noise sigma can be computed jointly for all clusters (default),
#' seperately for each cluster or be specified as input.
#' 
#' 
#' @param X a matrix with only continuous variables
#' @param segmentation a vector, segmentation for which likelihood is computed. Clusters
#'        numbers should be from range [1, numb.clusters]
#' @param max.dim an integer, maximum dimension of subspace. Number of principal components
#'        that span each subspace.
#' @param numb.clusters an integer, number of clusters
#' @param sigma a numeric, (default is NULL) value of sigma provided by the user
#' @param estimateJointly a boolean, (default value is TRUE) indicating if sigma should be estimated jointly for all clusters
#' @keywords internal
#' @return BIC value of BIC criterion
cluster.BIC <- function(X, segmentation, max.dim, numb.clusters, sigma=NULL, estimateJointly=TRUE){
  if(!is.matrix(X)){ # if X is one variable it is stored as vector
    X <- matrix(X, ncol=1)
  }
  D = dim(X)[1]
  p = dim(X)[2]
  if(is.null(sigma) & estimateJointly){
    RES.sigma <- sum(vapply(1:numb.clusters, function(k) {
      Xk = X[,segmentation==k, drop=F];
      if(dim(Xk)[2]>max.dim){ #length because it might be onedimensional
        svdSIGNAL= svd(Xk); 
        SIGNAL = matrix(svdSIGNAL$u[, 1:max.dim], ncol=max.dim) %*% 
          diag(svdSIGNAL$d[1:max.dim], nrow=max.dim) %*% 
          t(matrix(svdSIGNAL$v[, 1:max.dim], ncol=max.dim));
        return(sum((Xk - SIGNAL)^2))
      }
      return(0)
    }, 0.9))  
    degrees.freedom <- D*p-p-numb.clusters*D*max.dim-p*max.dim+numb.clusters*max.dim^2+numb.clusters*max.dim
    sigma <- sqrt(RES.sigma/degrees.freedom)
  }
  likelihoods <- rep(0, numb.clusters)  
  penalties <- rep(0, numb.clusters)  
  for(k in 1:numb.clusters){
    #one cluster
    Xk = X[,segmentation==k, drop=F]
    if(dim(Xk)[2]>max.dim){ #length because it might be onedimensional
      svdSIGNAL= svd(Xk)  
      SIGNAL = matrix(svdSIGNAL$u[, 1:max.dim], ncol=max.dim) %*% 
        diag(svdSIGNAL$d[1:max.dim], nrow=max.dim) %*% 
        t(matrix(svdSIGNAL$v[, 1:max.dim], ncol=max.dim))
      RESIDUAL = Xk - SIGNAL
      if(!estimateJointly & is.null(sigma)){ sigma = sqrt(sum(RESIDUAL^2)/((D-1)*(ncol(Xk)-1))) }
      likelihoods[k] <- sum(dnorm(as.matrix((RESIDUAL[,]), nrow=1), mean=0 , sd=sigma, log=T))
      mk <- ncol(Xk)
      penalties[k] <- 1/2*log(mk)*(max.dim*(D - max.dim +mk))
    }
    else{
      RESIDUAL = Xk - Xk
      likelihoods[k] <- sum(dnorm(as.matrix((RESIDUAL[,]), nrow=1), mean=0 , sd=sigma, log=T))
      mk <- max(ncol(Xk),1)
      penalties[k] <- 1/2*log(mk)*(max.dim*(D - max.dim +mk))
    }
  }
  BIC <- sum(likelihoods) - sum(penalties)
  return(BIC)
}


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


#' Plot mlcc.fit class object
#' 
#' @param x mlcc.fit class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
plot.mlcc.fit <- function(x,...){
  clusterNumbs <- lapply(x$all.fit, function(y) y$nClusters)
  BICVals <- lapply(x$all.fit, function(y) y$BIC)
  plot.default(clusterNumbs, BICVals, type="b", xaxt="n", ylab="BIC", xlab="Number of clusters")
  axis(side = 1, labels = clusterNumbs, at=clusterNumbs)
}

#' Print mlcc.fit class object
#' 
#' @param x mlcc.fit class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
print.mlcc.fit <- function(x,...){
  cat("$nClusters: ", x$nClusters, "\n")
  cat("$segmentation:\n")
  print(x$segmentation)
  cat("$BIC: ", x$BIC, "\n")
  cat("$subspacesDimensions:\n", unlist(x$subspacesDimensions), "\n")
}

#' Print mlcc.reps.fit class object
#' 
#' @param x mlcc.reps.fit class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
print.mlcc.reps.fit <- function(x,...){
  cat("$segmentation:\n")
  print(x$segmentation)
  cat("$BIC: ", x$BIC, "\n")
  cat("$basis:\n")
  cat(str(x$basis))
}