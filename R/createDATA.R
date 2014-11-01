#' Simulating subspace clustering data
#' 
#'  Generating data for simulation with a low-rank subspace structure: 
#'  variables are clustered and each cluster has a low-rank representation
#' 
#' @param n an integer, number of individuals
#' @param SNR a numeric, signal to noise ratio measured as variance of the variable, element of a subspace,
#'        to the variance of noise
#' @param K an integer, number of subspaces
#' @param numb.vars an integer, number of variables in each subspace
#' @param max.dim an integer, if equal.dims is TRUE then max.dim is dimension of each subspace.
#'        If equal.dims is FALSE then subspaces dimensions are drawn from uniform distribution on [1,max.dim]
#' @param equal.dims a boolean, if TRUE (value set by default) all clusters are of the same dimension
#' @export
#' @return A list consisting of:
#' \item{X}{matrix, generated data}
#' \item{signals}{matrix, data without noise}
#' \item{dims}{vector, dimensions of subspaces}
#' \item{s}{vector, true partiton of variables}
#' @examples
#' data <- data.simulation()
#' data2 <- data.simulation(n = 30, SNR = 2, K = 5, numb.vars = 20, max.dim = 3, equal.dims = FALSE)
data.simulation <- function(n = 100, SNR = 1, K = 10, numb.vars = 30, max.dim = 2, equal.dims = TRUE){
  sigma <- 1/SNR
  #subspaces dimensions depend on equal.dims value
  if(equal.dims)
    dims <- rep(max.dim,K)
  else
    dims <- sample(1:max.dim, K, replace=T)   
  
  X <- NULL
  Y <- NULL
  s <- NULL
  for (j in 1:K){
    Z <- qr.Q(qr(replicate(dims[j], rnorm(n, 0, 1))))
    coeff <- matrix(runif(dims[j]*numb.vars, 0.1, 1) * sign(runif(dims[j]*numb.vars, -1, 1)), nrow=dims[j])
    SIGNAL <- Z %*% coeff
    SIGNAL <- scale(SIGNAL)
    Y <- cbind(Y,SIGNAL)
    X <- cbind(X, SIGNAL + replicate(numb.vars, rnorm(n, 0, sigma)))
    s <- c(s, rep(j, numb.vars))
  }
  return(list(X = X,
              signals = Y,
              dims = dims,
              s = s))
}
