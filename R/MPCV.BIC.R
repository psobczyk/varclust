#' Subspace clustering based on multiple principal components
#'
#' Performs MPCC multiple times and chooses the best run based on residual variance
#'
#' @param X data
#' @param numb.Clusters clusters numbers to check
#' @param numb.runs number of runs of MLCC
#' @param stop.criterion how many changes in partitions triggers stopping the algorithm
#' @param max.iter maxium number of iteratations
#' @param max.dim maximum considered dimension of subspaces
#' @param method method to be used to determine best run. Possible values are "likelihood", "singular", "residual"
#' @param scale Should data be scaled?
#' @param numbCores Number of cores to be used
#' @param greedy Should number of clusters be estimated in a greedy way
#' @param estimateDimensions Should subspaces dimensions be estimated as well?
#' @export
#' @return a list consisting of
#' \item{segmentation}{of points to clusters}
#' \item{BIC}{Value of \code{\link{myBIC}} criterion}
MPCV.BIC <- function(X, numb.Clusters=1:10, numb.runs=20, stop.criterion=1, max.iter=20, max.dim=1, 
                     method=c("likelihood", "singular", "residual"), scale=T, numbCores=1, greedy=TRUE, estimateDimensions=T){
  if(scale){
    dane = scale(X)
  }
  method <- match.arg(method)
  n=nrow(X)
  p=ncol(X)
  results <- list()
  for(i in 1:length(numb.Clusters)){
    numb.clusters <- numb.Clusters[i]                                                                                                 
    MPCV.fit <- MPCV.reps(X=X, numb.Clusters=numb.clusters, numb.runs=numb.runs, max.dim=max.dim, method=method, scale=F, numbCores=numbCores)
    BIC.sum <- 0
    if(estimateDimensions){
      sigma <- NULL
      #SIGMA estimated jointly
      sigma <- getSigma(X, MPCV.fit$segmentation, max.dim, n, p, numb.clusters)
      #SIGMA estimated jointly 
      dimensions <- list() 
      for(k in 1:numb.clusters){
        temp <- 1
        for(d in 1:max.dim){ 
          temp[d] <- myBIC(X[,MPCV.fit$segmentation==k], rep(1, sum(MPCV.fit$segmentation==k)), d, 1, sigma=sigma)
        } 
        BIC.sum <- BIC.sum + max(temp[!is.nan(temp)]) 
        dimensions <- append(dimensions, which.max(temp)) 
      } 
      results[[i]] <- list(segmentation=MPCV.fit$segmentation, BIC=BIC.sum, subspacesDimensions=dimensions, nClusters=numb.clusters)
    }
    else{
      results[[i]] <- list(segmentation=MPCV.fit$segmentation, BIC=MPCV.fit$BIC, subspacesDimensions=NULL, nClusters=numb.clusters)
    }
    if(greedy & (i>2)){ 
      if( (results[[i]]$BIC < results[[i-2]]$BIC) & 
          (results[[i-1]]$BIC < results[[i-2]]$BIC)){
        break
      }
    }
  }
  BICs <- lapply(results, function(res) res$BIC)
  results[[which.max(BICs)]]
  
}
