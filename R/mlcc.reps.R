#' Multiple Latent Components Clustering - Subspace clustering assuming that the
#' number of clusters is known
#' 
#' For a fixed number of cluster function returns the best partition and basis 
#' for each subspace.
#' 
#' In more detail, an algorithm \code{\link{mlcc.kmeans}} is run a 
#' \emph{numb.runs} of times with random or custom initializations. The best 
#' partition is selected according to the BIC.
#' 
#' 
#' @param X A data frame or a matrix with only continuous variables.
#' @param numb.clusters An integer, number of cluster.
#' @param numb.runs An integer, number of runs of \code{\link{mlcc.kmeans} 
#'   algorithm} with random initialization.
#' @param stop.criterion An integer, if an iteration of 
#'   \code{\link{mlcc.kmeans}} algorithm makes less changes in partitions than 
#'   \code{stop.criterion}, \code{\link{mlcc.kmeans}} stops.
#' @param max.iter max.iter An integer, maximum number of iterations of the loop
#'   in \code{\link{mlcc.kmeans}} algorithm.
#' @param initial.segmentations A list of vectors, segmentations that user wants
#'   to be used as an initial segmentation in \code{\link{mlcc.kmeans}} 
#'   algorithm.
#' @param max.dim An integer, maximal dimension of subspaces.
#' @param scale A boolean, if TRUE (value set by default) then variables in 
#'   dataset are scaled to zero mean and unit variance.
#' @param numb.cores An integer, number of cores to be used, by default all 
#'   cores are used.
#' @param estimate.dimensions A boolean, if TRUE (value set by default) 
#'   subspaces dimensions are estimated.
#' @param flat.prior A boolean, if TRUE then, instead of a prior that takes into
#'   account number of models for a given number of clusters, flat prior is 
#'   used.
#' @param show.warnings A boolean, if set to TRUE all warnings are displayed, 
#'   default value is FALSE.
#' @param seed An integer, a seed for random number generator which allows 
#' making results reproductible
#' @export
#' @return A list consisting of \item{segmentation}{a vector containing the 
#'   partition of the variables} \item{BIC}{a numeric, value of the mBIC} 
#'   \item{basis}{a list of matrices, the factors for each of the subspaces}
#' @examples
#' \donttest{
#' sim.data <- data.simulation(n = 50, SNR = 1, K = 5, numb.vars = 50, max.dim = 3)
#' mlcc.reps(sim.data$X, numb.clusters = 5, numb.runs = 20, max.dim = 4)
#' }
mlcc.reps <- function(X, numb.clusters = 2, numb.runs = 30, stop.criterion = 1, max.iter = 30, 
  initial.segmentations = NULL, max.dim = 4, scale = TRUE, numb.cores = NULL, estimate.dimensions = TRUE, 
  flat.prior = FALSE, show.warnings = FALSE, seed = NULL) {
  if (is.data.frame(X)) {
    warning("X is not a matrix. Casting to matrix.")
    X <- as.matrix(X)
  }
  if (any(is.na(X))) {
    warning("Missing values are imputed by the mean of the variable")
    X[is.na(X)] = matrix(apply(X, 2, mean, na.rm = TRUE), ncol = ncol(X), nrow = nrow(X), 
      byrow = TRUE)[is.na(X)]
  }
  if (any(!sapply(X, is.numeric))) {
    auxi <- NULL
    for (j in 1:ncol(X)) if (!is.numeric(X[, j])) {
      auxi <- c(auxi, j)
    }
    stop(paste("\nThe following variables are not quantitative: ", auxi))
  }
  if (is.null(numb.cores)) {
    numb.cores <- max(1, detectCores() - 1)
  }
  cl <- makeCluster(numb.cores)
  registerDoParallel(cl)
  
  if (scale) {
    X <- scale(X)
  }
  
  i <- NULL
  BICs <- NULL
  segmentations <- NULL
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(initial.segmentations)) {
    segmentations <- foreach(i = (1:numb.runs)) %dorng% {
      MLCC.res <- mlcc.kmeans(X = X, number.clusters = numb.clusters, max.subspace.dim = max.dim, 
        max.iter = max.iter, estimate.dimensions = estimate.dimensions, show.warnings = show.warnings)
      current.segmentation <- MLCC.res$segmentation
      current.pcas <- MLCC.res$pcas
      list(current.segmentation, cluster.pca.BIC(X, current.segmentation, sapply(current.pcas, 
        ncol), numb.clusters, max.dim = max.dim, flat.prior = flat.prior), 
        current.pcas)
    }
    # running user specified clusters
  } else {
    segmentations <- foreach(i = (1:length(initial.segmentations))) %dorng% {
      MLCC.res <- mlcc.kmeans(X = X, number.clusters = numb.clusters, max.subspace.dim = max.dim, 
        max.iter = max.iter, initial.segmentation = initial.segmentations[[i]], 
        estimate.dimensions = estimate.dimensions, show.warnings = show.warnings)
      current.segmentation <- MLCC.res$segmentation
      current.pcas <- MLCC.res$pcas
      list(current.segmentation, cluster.pca.BIC(X, current.segmentation, sapply(current.pcas, 
        ncol), numb.clusters, max.dim = max.dim, flat.prior = flat.prior), 
        current.pcas)
    }
  }
  stopCluster(cl)
  BICs <- unlist(lapply(segmentations, function(x) x[2]))
  basis <- lapply(segmentations, function(x) x[3])
  segmentations <- lapply(segmentations, function(x) x[[1]])
  result <- list(segmentation = segmentations[[which.max(BICs)]], BIC = BICs[which.max(BICs)], 
    basis = basis[[which.max(BICs)]][[1]])
  class(result) <- "mlcc.reps.fit"
  return(result)
}
