#' Estimating number of relevant (non noise) Principal Components in PCA
#' 
#' For a given numeric data set, estimates number of PCs according to 
#' penalized likelihood criterion
#' 
#' 
#' @param X a data frame or a matrix contatining only continuous variables
#' @param ncp.min minimal number of principal components, for all the possible number of
#' PCs between ncp.min and ncp.max criterion is computed
#' @param ncp.max maximal number of principal components, if greater than dimensions of X,
#' min(ncol(X), nrow(X))-1 is used, for all the possible number of
#' PCs between ncp.min and ncp.max criterion is computed
#' @param scale a boolean, if TRUE (default value) then data is scaled before applying
#' criterion
#' @param verbose a boolean, if TRUE plot with BIC values for different
#'        numbers of components is produced (default value is FALSE)
#' @export
#' @return number of components
#' @examples
#' \dontrun{
#' library(MetabolAnalyze)
#' data(UrineSpectra)
#' estim.ncp(UrineSpectra[[1]], verbose=TRUE)}
estim.ncp <- function(X, ncp.min = 1, ncp.max = 10, scale = TRUE, verbose = FALSE){
  # preprocessing on X
  # number of components must be smaller than dimensions of X
  n <- nrow(X)
  p <- ncol(X)
  ncp.max <- min(ncp.max, min(n,p)-1)
  
  if(class(X) == "data.frame"){
    X <- as.matrix(X)
  }
  
  if(sum(sapply(X, is.numeric)) < p){
    stop("All the variables have to be numeric")
  }
  
  missing <- which(is.na(X))
  if(length(missing) !=  0){
    stop("There are missing values")
  }
  
  if(scale)
    X <- scale(X)
  
  # choosing the method to use accordingly to the dimensions of X
  method <- if(n>p){
    "Penalized likelihood, random factors model"
    } else "Penalized likelihood, random coefficients model"
  
  vals <- switch(method,
                 "Penalized likelihood, random coefficients model" = sapply(ncp.min:ncp.max, function(j) pca.new.BIC(X, j)),
                 "Penalized likelihood, random factors model" = sapply(ncp.min:ncp.max, function(j) pca.BIC(X, j)))
  if(verbose){
    caption <- paste0("Criterion:\n", method)
    plot(ncp.min:ncp.max, vals, xlab = "Number of components", ylab = "Criterion value",
         main = caption, type = "b")
    points(ncp.min-1+which.max(vals), max(vals), col = "red")
  }
  ncp.min-1+which.max(vals)
}
