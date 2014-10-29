#' Subspace clustering assuming that the number of clusters is known
#'
#' For a fixed number of cluster and fixed number of components per cluster, 
#' select the best partition according to the BIC.
#'
#' @param X a data frame or a matrix with only continuous variables
#' @param numb.clusters an integer, number of cluster
#' @param numb.runs an integer, number of runs of \code{\link{MPCV} algorithm} with random initialization
#' @param stop.criterion an integer, indicating how many changes in partitions triggers stopping the MLCC algorithm
#' @param max.iter an integer, maximum number of iterations of MLCC algorithm
#' @param initial.segmentations a list of vectors, segmentations that user wants to be used as an initial segmentation in MLCC algorithm
#' @param max.dim an integer, maximum dimension of subspaces to be considered
#' @param scale a boolean, if TRUE (value set by default) then data are scaled to unit variance
#' @param numb.cores an integer, number of cores to be used, by default all cores are used
#' @export
#' @return A list consisting of
#' \item{segmentation}{a vector containing the partition of the variables}
#' \item{BIC}{double, value of \code{\link{cluster.BIC}} criterion}
#' @examples
#' \donttest{
#' data <- dataSIMULATION(n=100, SNR=1, K=5, numbVars=30, max.dim=2)
#' MPCV.reps(data$X, numb.clusters=5, numb.runs=20)}
MPCV.reps <- function(X, numb.clusters=2, numb.runs=20, stop.criterion=1, max.iter=20, initial.segmentations=NULL,
                      max.dim=2, scale=T, numb.cores=NULL){
  if (is.data.frame(X)) {
    warnings("X is not a matrix. Casting to matrix.")
    X = as.matrix(X)
  }
  if(any(is.na(X))) {
    warnings("Missing values are imputed by the mean of the variable")
    X[is.na(X)] = matrix(apply(X, 2, mean, na.rm = TRUE), ncol = ncol(X), nrow = nrow(X), byrow = TRUE)[is.na(X)]
  }
  if (any(!sapply(X, is.numeric))) {
    auxi = NULL
    for (j in 1:ncol(X)) if (!is.numeric(X[, j])) 
      auxi = c(auxi, j)
    stop(paste("\nThe following variables are not quantitative: ", 
               auxi))
  }
  if (is.null(numb.cores)) {
    registerDoMC(max(1,detectCores()-1))
  }
  else{
    registerDoMC(numb.cores)
  }
  
  if(scale){
    X = scale(X)
  }
  else X=X
  
  BICs <- NULL 
  segmentations <- NULL
  segmentations <- foreach(i=(1:numb.runs)) %dopar% {
    MPCV.res <- MPCV(X=X, numberClusters=numb.clusters, maxSubspaceDim=max.dim, max.iter=max.iter)
    current.segmentation <- MPCV.res$segmentation
    current.pcas <- MPCV.res$pcas
    list(current.segmentation, cluster.BIC(X, current.segmentation, max.dim, numb.clusters))
  }
  i <- NULL
  segmentations2 <- foreach(i=(1:length(initial.segmentations))) %dopar% { #running user specified clusters
    MPCV.res <- MPCV(X=X, numberClusters=numb.clusters, maxSubspaceDim=max.dim, max.iter=max.iter, initial.segmentation=initial.segmentations[[i]])
    current.segmentation <- MPCV.res$segmentation
    current.pcas <- MPCV.res$pcas
    list(current.segmentation, cluster.BIC(X, current.segmentation, max.dim, numb.clusters))
  }
  segmentations <- append(segmentations, segmentations2)
  BICs <- unlist(lapply(segmentations, function(x) x[2]))
  segmentations <- lapply(segmentations, function(x) x[[1]])
  likelihood = return(list(segmentation = segmentations[[which.max(BICs)]], BIC = BICs[which.max(BICs)]))
}
