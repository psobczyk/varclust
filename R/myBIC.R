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
#' @param sigma a double, (default is NULL) value of sigma provided by the user
#' @param estimateJointly a boolean, (default value is TRUE) indicating if sigma should be estimated jointly for all clusters
#' @keywords internal
#' @return BIC value of BIC criterion
myBIC <- function(X, segmentation, max.dim, numb.clusters, sigma=NULL, estimateJointly=TRUE){
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


#' BIC for subspace clustering
#' 
#' Computes the value of BIC criterion for one cluster.
#' 
#' @param X data
#' @param max.dim maximum dimension of subspace
#' @param numb.clusters number of clusters
#' @param sigma (optional) pre-computed value of sigma
#' @keywords internal
#' @return BIC value of BIC criterion
myBIC_one_cluster <- function(X, max.dim, numb.clusters, sigma=NULL){
  if (!is.matrix(X)) {
    stop("In function myBIC_one_cluster. X is not a matrix.")
  }
  D = dim(X)[1]
  p = dim(X)[2]
  if(is.null(sigma)){
    RES.sigma <- 0
    Xk = X;
    if(length(unlist(Xk))>max.dim*D){ #length because it might be onedimensional
      svdSIGNAL= svd(Xk); 
      SIGNAL = matrix(svdSIGNAL$u[, 1:max.dim], ncol=max.dim) %*% 
        diag(svdSIGNAL$d[1:max.dim], nrow=max.dim) %*% 
        t(matrix(svdSIGNAL$v[, 1:max.dim], ncol=max.dim));
      RES.sigma <- (sum((Xk - SIGNAL)^2))
    }
    RES.sigma <- (0)
    degrees.freedom <- D*p-p-D*max.dim-p*max.dim+max.dim^2+max.dim
    sigma <- sqrt(RES.sigma/degrees.freedom)
  }
  
  Xk = X
  if(length(unlist(Xk))>=max.dim*D){ #length because it might be onedimensional
    svdSIGNAL= svd(Xk)  
    #print(sqrt(sum(svdSIGNAL$d[1:max.dim])))
    SIGNAL = matrix(svdSIGNAL$u[, 1:max.dim], ncol=max.dim) %*% 
      diag(svdSIGNAL$d[1:max.dim], nrow=max.dim) %*% 
      t(matrix(svdSIGNAL$v[, 1:max.dim], ncol=max.dim))
    RESIDUAL = Xk - SIGNAL
    #sigma = sqrt(sum(RESIDUAL^2)/((D-1)*(ncol(Xk)-1)))
    likelihood <- sum(dnorm(as.matrix((RESIDUAL[,]), nrow=1), mean=0 , sd=sigma, log=T))
  }
  #penalty on all clusters
  penalty = 1/2*log(p)*(1+ max.dim*(D + p - max.dim))
  BIC <- likelihood - penalty
  return(BIC)
}