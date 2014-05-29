#' Choose subspace closest to the given variable
#' 
#' The most similar subspace is choosen based on R^2
#' 
#' @param varaible
#' @param pcas orthogonal basis for different subspaces
#' @param numberClusters number of subspaces (clusters)
#' @return index number of subspace closest to variable
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