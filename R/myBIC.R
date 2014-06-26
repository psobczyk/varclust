#' BIC for subspace clustering
#' 
#' Computes the value of BIC criterion
#' 
#' @param X data
#' @param segmentation segmentation for which likelihood is computed
#' @param max.dim maximum dimension of subspace
#' @param numb.clusters number of clusters
#' @return BIC value of BIC criterion
myBIC <- function(X, segmentation, max.dim, numb.clusters){
  D = dim(X)[1]
  p = dim(X)[2]
  RES.sigma <- 0
  for(k in 1:numb.clusters){
    #one cluster
    Xk = X[,segmentation==k]
    if(length(unlist(Xk))>max.dim*D){ #length because it might be onedimensional
      svdSIGNAL= svd(Xk)  
      #print(sqrt(sum(svdSIGNAL$d[1:max.dim])))
      SIGNAL = matrix(svdSIGNAL$u[, 1:max.dim], ncol=max.dim) %*% 
        diag(svdSIGNAL$d[1:max.dim], nrow=max.dim) %*% 
        t(matrix(svdSIGNAL$v[, 1:max.dim], ncol=max.dim))
      RES.sigma = RES.sigma + sum((Xk - SIGNAL)^2)
    }
  }
  sigma <- sqrt(RES.sigma/((D-1)*(p-1)))
  likelihoods <- rep(0, numb.clusters)  
  for(k in 1:numb.clusters){
    #one cluster
    Xk = X[,segmentation==k]
    if(length(unlist(Xk))>=max.dim*D){ #length because it might be onedimensional
      svdSIGNAL= svd(Xk)  
      #print(sqrt(sum(svdSIGNAL$d[1:max.dim])))
      SIGNAL = matrix(svdSIGNAL$u[, 1:max.dim], ncol=max.dim) %*% 
               diag(svdSIGNAL$d[1:max.dim], nrow=max.dim) %*% 
               t(matrix(svdSIGNAL$v[, 1:max.dim], ncol=max.dim))
      RESIDUAL = Xk - SIGNAL
      #sigma = sqrt(sum(RESIDUAL^2)/((D-1)*(ncol(Xk)-1)))
      likelihoods[k] <- sum(dnorm(as.matrix((RESIDUAL[,]), nrow=1), mean=0 , sd=sigma, log=T))
    }
  }
  #penalty on all clusters
  penalty = 1/2*log(p)*(numb.clusters+ max.dim*(numb.clusters*D - numb.clusters*max.dim +p))
  BIC <- sum(likelihoods) - penalty
  return(BIC)
}