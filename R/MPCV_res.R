#' Subspace clustering based on multiple principal components
#'
#' Performs MPCC multiple times and chooses the best run based on residual variance
#'
#' @param X data
#' @param numb.Clusters number of clusters to be used
#' @param numb.runs number of runs of MLCC
#' @param stop.criterion how many changes in partitions triggers stopping the algorithm
#' @param max.iter maxium number of iteratations
#' @param initial.segmentation initial segmentation of variable to clusters
#' @param max.dim maximum considered dimension of subspaces
#' @param method method to be used to determine best run. Possible values are "likelihood", "singular", "residual"
#' @param scale Should data be scaled?
#' @param numbCores Number of cores to be used
#' @export
#' @return a list consisting of
#' \item{segmentation}{of points to clusters}
#' \item{BIC}{Value of \code{\link{myBIC}} criterion}
MPCV.reps <- function(X, numb.Clusters=2, numb.runs=20, stop.criterion=1, max.iter=20, initial.segmentation=NULL,
                      max.dim=1, method=c("likelihood", "singular", "residual"), scale=T, numbCores=1){
  registerDoMC(numbCores)
  
  method <- match.arg(method)
  if(scale){
    dane = scale(X)
  }
  else dane=X
  
  #goodness of fit
  Hs <- NULL
  Res <- NULL
  BICs <- NULL 
  
  segmentations <- NULL
  segmentations <- foreach(i=(1:numb.runs)) %dopar% {
    MPCV.res <- MPCV(X=dane, numberClusters=numb.Clusters, maxSubspaceDim=max.dim, max.iter=max.iter, initial.segmentation=initial.segmentation)
    current.segmentation <- MPCV.res$segmentation
    current.pcas <- MPCV.res$pcas
    
    if(method=="singular"){
      H <- 0 #explained variance
      for (j in 1:numb.Clusters){
        if(is.matrix(dane[,current.segmentation==j]) && dim(dane[,current.segmentation==j])[2]>max.dim){
          a <- summary(prcomp(dane[,current.segmentation==j]))
          H <- H + a$importance[3,max.dim]
        }
        else{
          H <- H+1 #all variance in that cluster is explained
        }
      }
      list(current.segmentation, myBIC(dane, current.segmentation, max.dim, numb.Clusters), H)
    }
    
    if(method=="residual"){
      R <- 0 #residual variance
      for (j in 1:numb.Clusters){
        if(is.matrix(dane[,current.segmentation==j]) && dim(dane[,current.segmentation==j])[2]>max.dim){
          a <- summary(prcomp(dane[,current.segmentation==j]))
          R <- R + sum(lm(dane[,current.segmentation==j] ~ current.pcas[[j]])$residuals^2)
        } #otherwise we have a perfect fit
      }
      list(current.segmentation, myBIC(dane, current.segmentation, max.dim, numb.Clusters), R)
    }
    else{
      list(current.segmentation, myBIC(dane, current.segmentation, max.dim, numb.Clusters))
    }
  }
  BICs <- unlist(lapply(segmentations, function(x) x[2]))
  keep <- segmentations
  segmentations <- lapply(segmentations, function(x) x[[1]])
  switch(method,
         singular   = return(list(segmentation = segmentations[[which.max(Hs)]],   BIC = BICs[which.max(unlist(lapply(keep, function(x) x[3])))])),
         residual   = return(list(segmentation = segmentations[[which.min(Res)]],  BIC = BICs[which.min(unlist(lapply(keep, function(x) x[3])))])),
         likelihood = return(list(segmentation = segmentations[[which.max(BICs)]], BIC = BICs[which.max(BICs)])))
}
