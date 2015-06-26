#' Estimating number of important Principal Components in PCA
#' 
#' For a given data set, estimates numberof PCs according to 
#' one of four different penalized likelihood criterions
#' Methods description:
#' \itemize{
#'  \item penalizedLikelihood - standard penalized likelihood
#'  \item homogenousPenalizedLikelihood - penalized likelihood with singular values assumed equal
#' and coefficients with non-one variance
#'  \item minkaBIC - BIC criterion dervied by Minka in Automatic choice of dimensionality for PCA
#'  \item laplaceEvidence - prior dependent criterion dervied by Minka in Automatic choice of dimensionality for PCA
#' }
#' 
#' @param X a data frame or a matrix with only continuous variables
#' @param ncp.min minimal number of principal components
#' @param ncp.max maximal number of principal components
#' @param method a crtierion to be used
#' @param scale a boolean, if TRUE (default value) then data is scaled
#' @param verbose a boolean, if TRUE plot with BIC values for different
#'        numbers of components is produced 
#' @export
#' @references P. Sobczyk, M. Bogdan, J. Josse, 
#' Bayesian dimensionality reduction with PCA using partially integrated penalized likelihood
#' @return Number of components
#' 
estim.ncp <- function(X, ncp.min = 1, ncp.max = 10, method = c("penalizedLikelihood", "homogenousPenalizedLikelihood", 
                                                               "minkaBIC", "laplaceEvidence"), 
                      scale = TRUE, verbose = FALSE){
  method <- match.arg(method)
  if(scale)
    X <- scale(X)
  vals <- switch(method,
                 penalizedLikelihood = sapply(ncp.min:ncp.max, function(j) pca.new.BIC(X, j)),
                 homogenousPenalizedLikelihood = sapply(ncp.min:ncp.max, function(j) rajan.BIC(X, j)),
                 minkaBIC = sapply(ncp.min:ncp.max, function(j) pca.BIC(X, j)),
                 laplaceEvidence = sapply(ncp.min:ncp.max, function(j) pca.Laplace(X, j)))
  if(verbose){
    caption <- paste0("Criterion: ", method)
    plot(ncp.min:ncp.max, vals, xlab = "Number of components", ylab = "Criterion value",
         main = caption, type = "b")
    points(ncp.min-1+which.max(vals), max(vals), col = "red")
  }
  ncp.min-1+which.max(vals)
}
