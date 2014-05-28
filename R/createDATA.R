#' Piotr Sobczyk IMiI
#' Generating random data for subspace clustering simulation
#' Result is matrix X

dataSIMULATION <- function(n = 100, SNR, K = 20, numbPoints = 50, max_dim = 1){
  #draw dimensions of subspaces
  sigma = 1/(SNR*sqrt(n*K*numbPoints))
  dims = rep(1, K) #sample(1:max_dim, K, replace=T) 
  
  X = NULL
  Y = NULL
  s = NULL
  for (j in 1:K){
    SIGNAL = replicate(numbPoints, rnorm(n, 0, 1))
    SIGNAL = scale(SIGNAL, scale = FALSE)
    svdSIGNAL= svd(SIGNAL)  
    SIGNAL = matrix(svdSIGNAL$u[, 1:dims[j]], ncol=dims[j]) %*% diag(svdSIGNAL$d[1:dims[j]], nrow=dims[j]) %*% 
      t(matrix(svdSIGNAL$v[, 1:dims[j]], ncol=dims[j])) / sqrt(sum(svdSIGNAL$d[1:dims[j]]^2))
    X = cbind(X, SIGNAL + sigma*replicate(numbPoints, rnorm(n, 0, 1)))
    Y = cbind(Y, SIGNAL)
    s = c(s, rep(j, numbPoints))
  }
  return(list(X = X,
              signals = Y,
              dims = dims,
              s = s))
}
