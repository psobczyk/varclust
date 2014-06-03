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
  likelihoods <- NULL  
  for(k in 1:numb.clusters){
    #one cluster
    Xk = X[,segmentation==k]
    if(length(Xk)>=max.dim*D){ #length because it might be onedimensional
      svdSIGNAL= svd(Xk)  
      SIGNAL = matrix(svdSIGNAL$u[, 1:max.dim], ncol=max.dim) %*% diag(svdSIGNAL$d[1:max.dim], nrow=max.dim) %*% 
        t(matrix(svdSIGNAL$v[, 1:max.dim], ncol=max.dim)) #/ sqrt(sum(svdSIGNAL$d[1:max.dim]^2))
      RESIDUAL = Xk - SIGNAL
      sigma = sqrt(sum(RESIDUAL^2)/(D*ncol(Xk)))
      likelihoods[k] <- sum(dnorm(as.matrix((RESIDUAL[,]), nrow=1), mean=0 , sd=sigma, log=T))
    }
  }
  #penalty on all clusters
  penalty = log(D)/2*max.dim*(numb.clusters*D + numb.clusters*max.dim +p)
  BIC <- sum(likelihoods) - penalty
  return(BIC)
}