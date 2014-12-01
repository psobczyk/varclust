#' Subspace clustering assuming that the number of clusters is known
#'
#' For a fixed number of cluster and fixed number of components per cluster
#' function returns the best partition and basis for each subspace.
#' 
#' In more detail, an algorithm \code{\link{mlcc.kmeans}} 
#' is run a \emph{numb.runs} of times with random initializations. 
#' The best partition is selected according to the BIC.
#' 
#'
#' @param X a data frame or a matrix with only continuous variables
#' @param numb.clusters an integer, number of cluster
#' @param numb.runs an integer, number of runs of \code{\link{mlcc.kmeans} algorithm} with random initialization
#' @param stop.criterion an integer, if an iteration of \code{\link{mlcc.kmeans}} algorithm 
#'        makes less changes in partitions than \code{stop.criterion}, 
#'        \code{\link{mlcc.kmeans}} stops.
#' @param max.iter an integer, maximum number of iterations of \code{\link{mlcc.kmeans}} algorithm
#' @param initial.segmentations a list of vectors, segmentations that user wants to be 
#'        used as an initial segmentation in \code{\link{mlcc.kmeans}} algorithm
#' @param max.dim an integer, dimension of subspaces (all are assumed to be equal)
#' @param scale a boolean, if TRUE (value set by default) then variables in 
#'        dataset are scaled to zero mean and unit variance
#' @param numb.cores an integer, number of cores to be used, by default all cores are used
#' @param estimate.dimensions a boolean, if TRUE (value set by default) subspaces dimensions are estimated
#' @export
#' @return A list consisting of
#' \item{segmentation}{a vector containing the partition of the variables}
#' \item{BIC}{a numeric, value of \code{\link{cluster.BIC}} criterion}
#' \item{basis}{a list of matrices, the basis vectors for subspaces}
#' @examples
#' \donttest{
#' sim.data <- data.simulation(n = 100, SNR = 1, K = 5, numb.vars = 30, max.dim = 2)
#' mlcc.reps(sim.data$X, numb.clusters = 5, numb.runs = 20, max.dim = 4)
#' }
mlcc.reps <- function(X, numb.clusters = 2, numb.runs = 20, stop.criterion = 1, max.iter = 20, 
                      initial.segmentations = NULL, max.dim = 4, scale = TRUE, numb.cores = NULL,
                      estimate.dimensions = TRUE){
  if (is.data.frame(X)) {
    warning("X is not a matrix. Casting to matrix.")
    X = as.matrix(X)
  }
  if (any(is.na(X))) {
    warning("Missing values are imputed by the mean of the variable")
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
    X <- scale(X)
  }
  else X <- X
  
  BICs <- NULL 
  segmentations <- NULL
  segmentations <- foreach(i=(1:numb.runs)) %dopar% {
    MPCV.res <- mlcc.kmeans(X=X, number.clusters=numb.clusters, max.subspace.dim=max.dim, max.iter=max.iter, 
                            estimate.dimensions = estimate.dimensions)
    current.segmentation <- MPCV.res$segmentation
    current.pcas <- MPCV.res$pcas
    list(current.segmentation, 
         adjusted.cluster.BIC(X, current.segmentation, 
                              sapply(current.pcas, ncol), numb.clusters), 
         current.pcas)
  }
  i <- NULL
  #running user specified clusters
  segmentations2 <- foreach(i=(1:length(initial.segmentations))) %dopar% { 
    MPCV.res <- mlcc.kmeans(X = X, number.clusters = numb.clusters, 
                            max.subspace.dim = max.dim, max.iter = max.iter, 
                            initial.segmentation = initial.segmentations[[i]],
                            estimate.dimensions = estimate.dimensions)
    current.segmentation <- MPCV.res$segmentation
    current.pcas <- MPCV.res$pcas
    list(current.segmentation, 
         adjusted.cluster.BIC(X, current.segmentation, 
                              sapply(current.pcas, ncol), numb.clusters), 
         current.pcas)
  }
  segmentations <- append(segmentations, segmentations2)
  BICs <- unlist(lapply(segmentations, function(x) x[2]))
  basis <- lapply(segmentations, function(x) x[3])
  segmentations <- lapply(segmentations, function(x) x[[1]])
  result <- list(
              segmentation = segmentations[[which.max(BICs)]], 
              BIC = BICs[which.max(BICs)],
              basis = basis[[which.max(BICs)]][[1]])
  class(result) <- "mlcc.reps.fit"
  return(result)
}
