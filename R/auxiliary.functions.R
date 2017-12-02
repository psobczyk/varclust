#' Penalized likelihood for PCA
#' 
#' Computes the value of BIC-like criterion for given data set and 
#' number of factors. Assumes that number of variables is large
#' compared to number of observations
#' 
#' @param X a matrix with only continuous variables
#' @param k number of principal components fitted
#' @keywords internal
#' @return BIC value of BIC criterion
pca.new.BIC <- function(X, k){
  d <- dim(X)[1]
  N <- dim(X)[2]
  m <- d*k - k*(k+1)/2
  
  ## TO DO: replace eigen with some more robust eigenvalues computation
  lambda <- eigen(cov(t(X)), only.values = TRUE)$values
  if(any(lambda < 0)){
    #warning("In function pca.new.BIC: some of the eigenvalues were negative due to numerical errors - rounding them to 0")
    lambda[lambda < 0] = 0
  }
  v <- sum(lambda[(k+1):d])/(d-k) 
  
  t0 <- -N*d/2*log(2*pi)
  t1 <- -N/2*sum(log(lambda[1:k]))
  t2 <- -N*(d-k)/2*log(v)
  t3 <- -N*d/2
  pen <- -(m+d+k+1)/2*log(N)
  t0+t1+t2+t3+pen
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
#' @param dims a vector of integers, dimensions of subspaces. Number of principal components
#'        that span each subspace.
#' @param numb.clusters an integer, number of clusters
#' @param max.dim an integer, upper bound for allowed dimension of subspace
#' @param flat.prior boolean, if TRUE (default is FALSE) then flat prior on models is used
#' @keywords internal
#' @return BIC value of BIC criterion
cluster.pca.BIC <- function(X, segmentation, dims, numb.clusters, max.dim, flat.prior = FALSE){
  if(!is.matrix(X)){ # if X is one variable it is stored as vector
    X <- matrix(X, ncol=1)
  }
  D = dim(X)[1]
  p = dim(X)[2]
  
  formula <- rep(0, numb.clusters)   
  for(k in 1:numb.clusters){
    #one cluster
    dim1 = dims[k]
    Xk = X[,segmentation==k, drop=F]
    ## TO DO: use R pesel package
    if(dim(Xk)[2] > dim1){
      if(dim(Xk)[2] > D){
        formula[k] <- pca.new.BIC(Xk, dim1)
      } else{
        formula[k] <- pca.new.BIC(t(Xk), dim1)
      }
    } else{
      ## TO DO: reconsider this!!! Possible undesired behaviour
      formula[k] <- - Inf
    }
  }
  #apriori
  apriori.segmentations <- -p*log(numb.clusters)
  apriori.dimensions <- - log(max.dim)*numb.clusters
  BIC <- sum(formula)
  if (!flat.prior){
    BIC <- BIC + apriori.segmentations + apriori.dimensions
  }
  return(BIC)
}


#' Assigns variables to subspaces (according to BIC)
#'
#' Selects subspace closest to given variable (according to BIC)
#'
#' @param variable variable variable to be assigned
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
    sigma.hat <- sqrt(sum(res^2)/n)
    if (sigma.hat < 1e-15){
      #TO DO: add parametr verbose so we get rid of warning every time we intialize
      #warning("In function choose.cluster.BIC: estimated value of noise in cluster is <1e-15. It might corrupt the result.")
    }
    loglik <- sum(dnorm(res, 0, sigma.hat, log=T))
    BICs[i] <- loglik - nparams*log(n)/2
  }
  which.max(BICs)
}

#' Calculates pseudo-distance from given variable to one dimensional clusters by BIC (initialization of kmeans++)
#'
#' @param variable variable for which the pseudo-distance is calculated 
#' @param pcas orthogonal basis for different subspaces
#' @param numberClusters number of subspaces (clusters)
#' @keywords internal
#' @return minimal distance to a subspace
calculate.distance.kmeanspp <- function(variable, pcas, numberClusters){
  dists <- NULL
  for(i in 1:numberClusters){
    if(ncol(pcas[[i]]) != 1){
      stop("For one dimensional clusters only")
    }
    n <- length(variable)
    res <- fastLmPure(pcas[[i]], variable, method = 0L)$residuals
    dists[i] <- sum(res^2)/n
  }
  min(dists)
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