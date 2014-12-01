#' Subspace clustering with automatic estimation of number of clusters
#' and their dimension
#'
#' Estimate the number of clusters according to the BIC. Basic k-means based 
#' Multiple Latent Components Clustering (MLCC) algorithm (\code{\link{mlcc.kmeans}}) is run a 
#' given number of times (\emph{numb.runs}) for each number of clusters in \emph{numb.clusters}.
#' The best partition is choosen with BIC (see \code{\link{mlcc.reps} function})
#'
#' @param X a data frame or a matrix with only continuous variables
#' @param numb.clusters a vector, numbers of clusters to be checked
#' @param numb.runs an integer, number of runs of \code{\link{mlcc.kmeans}}
#' @param stop.criterion an integer, if an iteration of \code{\link{mlcc.kmeans}} algorithm 
#'        makes less changes in partitions than \code{stop.criterion}, 
#'        \code{\link{mlcc.kmeans}} stops.
#' @param max.iter an integer, maximum number of iterations of \code{\link{mlcc.kmeans}} algorithm
#' @param max.dim an integer, if estimate.dimensions is FALSE then max.dim is dimension of each subspace.
#'        If estimate.dimensions is TRUE then subspaces dimensions are estimated from the range [1, max.dim]
#' @param scale a boolean, if TRUE (value set by default) then variables in 
#'        dataset are scaled to zero mean and unit variance
#' @param numb.cores an integer, number of cores to be used, by default all cores are used
#' @param greedy a boolean, if TRUE (value set by default) the clusters are estimated in a greedy way
#' @param estimate.dimensions a boolean, if TRUE (value set by default) subspaces dimensions are estimated
#' @param verbose a boolean, if TRUE plot with BIC values for different
#'        numbers of clusters is produced and values of BIC, computed
#'        for every number of clusters and subspaces dimensions, are printed
#'        (value set by default is FALSE)
#' @export
#' @return An object of class mlcc.fit consisting of
#' \item{segmentation}{a vector containing the partition of the variables}
#' \item{BIC}{numeric, value of \code{\link{cluster.BIC}} criterion}
#' \item{subspacesDimensions}{a list containing dimensions of the subspaces}
#' \item{nClusters}{an integer, estimated number of clusters}
#' \item{factors}{a list of matrices, basis for each subspace}
#' \item{all.fit}{a list of segmentation, BIC, subspaces dimension for all numbers of clusters considered for an estimated subspace dimensions}
#' \item{all.fit.dims}{a list of lists of segmentation, BIC, subspaces dimension for all numbers of clusters and subspaces dimensions considered}
#' @examples
#' \donttest{
#' sim.data <- data.simulation(n = 100, SNR = 1, K = 5, numb.vars = 30, max.dim = 2)
#' mlcc.bic(sim.data$X, numb.clusters = 1:10, numb.runs = 20, verbose=TRUE)
#' }
mlcc.bic <- function(X, numb.clusters = 1:10, numb.runs = 20, stop.criterion = 1, max.iter = 20, max.dim = 4, 
                    scale = TRUE, numb.cores = NULL, greedy = TRUE, estimate.dimensions = TRUE, verbose = FALSE){
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
  if (scale){
    X = scale(X)
  }
  n <- nrow(X)
  p <- ncol(X)
  greedy.stop <- max(numb.clusters)
  results <- list()
  if (verbose)
    cat("Number of clusters \t BIC \n")
  
  for(i in 1:length(numb.clusters)){
    number.clusters <- numb.clusters[i]
    MPCV.fit <- mlcc.reps(X = X, numb.clusters = number.clusters, 
                          numb.runs = numb.runs, max.dim = max.dim, scale = F,
                          numb.cores = numb.cores, 
                          estimate.dimensions = estimate.dimensions)
    
    results[[i]] <- list(segmentation = MPCV.fit$segmentation, 
                         BIC = MPCV.fit$BIC, 
                         subspacesDimensions = lapply(MPCV.fit$basis, ncol), 
                         nClusters = number.clusters, factors = MPCV.fit$basis)
    if (greedy & (i>2)) {
      if ( (results[[i]]$BIC < results[[i-1]]$BIC) &
             (results[[i-2]]$BIC < results[[i-1]]$BIC) ){
        greedy.stop <- i
        if (verbose) cat(paste("       ", number.clusters, "        ", 
                               formatC(results[[i]]$BIC, digits = 
                                         ceiling(log(abs(results[[i]]$BIC),10))), "\n"))
        break
      }
    }
    if (verbose) cat(paste("       ", number.clusters, "        ",
                           formatC(results[[i]]$BIC, digits = 
                                     ceiling(log(abs(results[[i]]$BIC),10))), "\n"))
  }
  BICs <- lapply(results, function(res) res$BIC)
  if (verbose) {
    plot(numb.clusters[1:greedy.stop], BICs, type="b", xaxt="n", ylab="BIC", 
         xlab="Number of clusters")
    axis(side = 1, labels = numb.clusters[1:greedy.stop], at=numb.clusters[1:greedy.stop])
  }
  result <- results[[which.max(BICs)]]
  result$factors <- lapply(1:result$nClusters, function(i) {
    d <- ncol(result$factors[[i]]); 
    colnames(result$factors[[i]]) <- paste(i, 1:d); 
    result$factors[[i]] } )
  result$all.fit <- results
  class(result) <- "mlcc.fit"
  return(result)
}
