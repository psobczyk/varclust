#' Subspace clustering with automatic estimation of number of clusters
#' and their dimension.
#'
#' Estimate the number of clusters according to the BIC. Basic k-means based 
#' Multiple Latent Components Clustering (MLCC) algorithm (\code{\link{mlcc.kmeans}}) is run a 
#' given number of times (\emph{numb.runs}) for each number of clusters in numb.Clusters.
#' The best partition is choosen with BIC (see \code{\link{mlcc.reps} function})
#'
#' @param X a data frame or a matrix with only continuous variables
#' @param numb.clusters a vector, numbers of clusters to be checked
#' @param numb.runs an integer, number of runs of \code{\link{mlcc.kmeans}}
#' @param stop.criterion an integer, indicating how many changes in partitions triggers stopping the \code{\link{mlcc.kmeans} algorithm}
#' @param max.iter an integer, maximum number of iterations of \code{\link{mlcc.kmeans} algorithm}
#' @param max.dim an integer, maximum dimension of subspaces to be considered
#' @param scale a boolean, if TRUE (value set by default) then data are scaled to unit variance
#' @param numb.cores an integer, number of cores to be used, by default all cores are used
#' @param greedy a boolean, if TRUE (value set by default) the clusters are estimated in a greedy way
#' @param estimate.dimensions a boolean, if TRUE (value set by default) subspaces dimensions are estimated
#' @param graphical.output a boolean, if TRUE (value set by default) plot with BIC values for different
#'        numbers of clusters is produced
#' @export
#' @return An object of class mlcc.fit consisting of
#' \item{segmentation}{a vector containing the partition of the variables}
#' \item{BIC}{double, value of \code{\link{cluster.BIC}} criterion}
#' \item{subspacesDimensions}{a list containing dimensions of the subspaces}
#' \item{nClusters}{an integer, estimated number of clusters}
#' \item{all.fit}{a list, segmentation, BIC, subspaces dimension for all numbers of clusters considered}
#' @examples
#' \donttest{
#' data <- data.simulation(n=100, SNR=1, K=5, numb.vars=30, max.dim=2)
#' mlcc.bic(data$X, numb.clusters=1:10, numb.runs=20)
#' }
mlcc.bic <- function(X, numb.clusters=1:10, numb.runs=20, stop.criterion=1, max.iter=20, max.dim=4, 
                    scale=TRUE, numb.cores=NULL, greedy=TRUE, estimate.dimensions=TRUE, graphical.output=TRUE){
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
  if(scale){
    X = scale(X)
  }
  n=nrow(X)
  p=ncol(X)
  greedy.stop <- max(numb.clusters)
  results <- list()
  cat(paste("Number of clusters \t BIC \n"))
  for(i in 1:length(numb.clusters)){
    number.clusters <- numb.clusters[i]                                                                                                 
    MPCV.fit <- mlcc.reps(X=X, numb.clusters=number.clusters, numb.runs=numb.runs, max.dim=max.dim, scale=F, numb.cores=numb.cores)
    BIC.sum <- 0
    if(estimate.dimensions){
      sigma <- NULL
      #SIGMA estimated jointly
      sigma <- getSigma(X, MPCV.fit$segmentation, max.dim, n, p, number.clusters)
      #SIGMA estimated jointly 
      dimensions <- list() 
      for(k in 1:number.clusters){
        temp <- 1
        for(d in 1:max.dim){ 
          temp[d] <- cluster.BIC(X[,MPCV.fit$segmentation==k, drop=F], rep(1, sum(MPCV.fit$segmentation==k)), d, 1, sigma=sigma)
        } 
        BIC.sum <- BIC.sum + max(temp[!is.nan(temp)]) 
        dimensions <- append(dimensions, which.max(temp)) 
      } 
      results[[i]] <- list(segmentation=MPCV.fit$segmentation, BIC=BIC.sum, subspacesDimensions=dimensions, nClusters=number.clusters)
    }
    else{
      results[[i]] <- list(segmentation=MPCV.fit$segmentation, BIC=MPCV.fit$BIC, subspacesDimensions=NULL, nClusters=number.clusters)
    }
    if(greedy & (i>2)){ 
      if ( (results[[i]]$BIC < results[[i-2]]$BIC) & 
           (results[[i-1]]$BIC < results[[i-2]]$BIC) ){
        greedy.stop <- i
        cat(paste("        ", number.clusters, "        ", formatC(results[[i]]$BIC, digits=ceiling(log(abs(results[[i]]$BIC),10))), "\n"))
        break
      }
    }
    cat(paste("        ", number.clusters, "        ", formatC(results[[i]]$BIC, digits=ceiling(log(abs(results[[i]]$BIC),10))), "\n"))
  }
  BICs <- lapply(results, function(res) res$BIC)
  if (graphical.output) {
    plot(numb.clusters[1:greedy.stop], BICs, type="b", xaxt="n", ylab="BIC", xlab="Number of clusters")
    axis(side = 1, labels = numb.clusters[1:greedy.stop], at=numb.clusters[1:greedy.stop])
  }

  result <- results[[which.max(BICs)]]
  result$all.fit <- results
  class(result) <- "mlcc.fit"
  return(result)
}
