require(RcppEigen)

choose.cluster <- function(variable, pcas, numberClusters){
  rSquare <- NULL
  v1 = var(variable)
  for(i in 1:numberClusters){
    v2 <- var(fastLmPure(pcas[[i]], variable, method = 0L)$residuals)
    p <- ncol(pcas[[i]])
    n <- length(variable)
    #rSquare[i] <- 1 - ( 1- (v1-v2)/v1) *(n-1)/(n-p-1)
    rSquare[i] <- (v1-v2)/v1
  }
  which.max(rSquare)
}


#BIC for subspace clustering
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
        t(matrix(svdSIGNAL$v[, 1:max.dim], ncol=max.dim)) / sqrt(sum(svdSIGNAL$d[1:max.dim]^2))
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