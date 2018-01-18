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
    dimk = dims[k]
    Xk = X[,segmentation==k, drop=F]
    if(dim(Xk)[2] > dimk && dim(Xk)[1] > dimk){
      formula[k] <- pesel(X = Xk, npc.min = dimk, npc.max = dimk, scale = FALSE, method = "heterogenous")$vals[1]
    } else{
      #if after reassignment of variables to clusters there are less variables in it than the calculated dimensionality of the cluster
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

#' Selects subspace closest to a given variable (according to BIC)
#'
#' @param variable variable variable to be assigned
#' @param pcas orthogonal basis for different subspaces
#' @param numberClusters number of subspaces (clusters)
#' @param show.warnings a boolean - if set to TRUE all warnings are displayed, default value is FALSE
#' @keywords internal
#' @return index number of subspace closest to variable
choose.cluster.BIC <- function(variable, pcas, numberClusters, show.warnings = FALSE){
  BICs <- NULL
  for(i in 1:numberClusters){
    nparams <- ncol(pcas[[i]])
    n <- length(variable)
    res <- fastLmPure(pcas[[i]], variable, method = 0L)$residuals
    sigma.hat <- sqrt(sum(res^2)/n)
    if (sigma.hat < 1e-15 && show.warnings){
      warning("In function choose.cluster.BIC: estimated value of noise in cluster is <1e-15. It might corrupt the result.")
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