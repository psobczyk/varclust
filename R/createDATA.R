#' Simulating subspace clustering data
#' 
#' Generating random data for subspace clustering simulation
#' 
#' @param n an integer, number of individuals
#' @param SNR a double, signal to noise ratio measured as variance of variable to variance of noise
#' @param K an integer, number of subspaces
#' @param numbVars an integer, number of variables in each subspace
#' @param max.dim an integer, maximum dimension of subspace
#' @param equalDims a boolean, if TRUE (value set by default) all clusters are of the same dimension
#' @export
#' @return A list consisting of:
#' \item{X}{matrix, generated data}
#' \item{signals}{matrix, data without noise}
#' \item{dims}{vector, dimensions of subspaces}
#' \item{s}{vector, true partiton of variables}
dataSIMULATION <- function(n = 100, SNR=1, K = 10, numbVars = 30, max.dim = 2, equalDims=T){
  #draw dimensions of subspaces
  sigma = 1/SNR
  if(equalDims)
    dims = rep(max_dim,K)
  else
    dims = sample(1:max_dim, K, replace=T)   
  
  X = NULL
  Y = NULL
  s = NULL
  for (j in 1:K){
    Z <- qr.Q(qr(matrix(mvrnorm(dims[j],rep(0,n),diag(rep(1,n))), ncol=dims[j])))
    coeff <- matrix(runif(dims[j]*numbPoints, 0.1, 1) * sign(runif(dims[j]*numbPoints, -1, 1)), nrow=dims[j])
    SIGNAL <- Z %*% coeff
    SIGNAL <- scale(SIGNAL)
    Y = cbind(Y,SIGNAL)
    X = cbind(X, SIGNAL + replicate(numbPoints, rnorm(n, 0, sigma)))
    s = c(s, rep(j, numbPoints))
  }
  return(list(X = X,
              signals = Y,
              dims = dims,
              s = s))
}
