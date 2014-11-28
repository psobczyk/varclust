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
#' @param dims a vector of integers, dimensions of subspaces. Number of principal components
#'        that span each subspace.
#' @param numb.clusters an integer, number of clusters
#' @param adjustment a numeric, percentage of BIC penalty applied
#' @param sigma a numeric, (default is NULL) value of sigma provided by the user
#' @param estimateJointly a boolean, (default value is TRUE) indicating if sigma should be estimated jointly for all clusters
#' @keywords internal
#' @return BIC value of BIC criterion
adjusted.cluster.BIC <- function(X, segmentation, dims, numb.clusters, adjustment = 1, sigma=NULL, estimateJointly=TRUE){
  if(!is.matrix(X)){ # if X is one variable it is stored as vector
    X <- matrix(X, ncol=1)
  }
  D = dim(X)[1]
  p = dim(X)[2]
  if(is.null(sigma) & estimateJointly){
    RES.sigma <- sum(vapply(1:numb.clusters, function(k) {
      max.dim = dims[k]
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
    max.dim = max(dims)
    degrees.freedom <- D*p-p-numb.clusters*D*max.dim-p*max.dim+numb.clusters*max.dim^2+numb.clusters*max.dim
    sigma <- sqrt(RES.sigma/degrees.freedom)
  }
  likelihoods <- rep(0, numb.clusters)  
  penalties <- rep(0, numb.clusters)  
  for(k in 1:numb.clusters){
    #one cluster
    max.dim = dims[k]
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
  BIC <- sum(likelihoods) - adjustment*sum(penalties)
  return(BIC)
}


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
  if (!is.matrix(X)){ # if X is one variable it is stored as vector
    X <- matrix(X, ncol=1)
  }
  D = dim(X)[1]
  p = dim(X)[2]
  if (is.null(sigma) & estimateJointly){
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
    if (dim(Xk)[2]>max.dim){ #length because it might be onedimensional
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

#' Selects subspace closest to given variable (according to BIC)
#'
#' The most similar subspace is choosen based on BIC criterion
#'
#' @param variable variable to be assigned
#' @param pcas orthogonal basis for different subspaces
#' @param numberClusters number of subspaces (clusters)
#' @keywords internal
#' @return index number of subspace closest to variable
choose.cluster.BIC <- function(variable, pcas, numberClusters){
  BICs <- NULL
  for(i in 1:numberClusters){
    nparams <- ncol(pcas[[i]])
    n <- length(variable)
    res <- fastLmPure(pcas[[i]], variable, method = 0L)$residuals
    sigma.hat <- sqrt(sum(res^2)/100)
    if (sigma.hat < 1e-15)
      warning("In function choose.cluster.BIC: estimated value of noise in cluster is <1e-16. It might corrupt the result.")
    loglik <- sum(dnorm(res, 0, sigma.hat, log=T))
    BICs[i] <- loglik - nparams*log(n)/2
  }
  which.max(BICs)
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