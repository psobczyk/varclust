#' Subspace clustering based on multiple principal components
#'
#' Performs MPCC multiple times and chooses the best run based on residual variance
#'
#' @param X data
#' @param numberClusters number of clusters to be used
#' @param numb.runs number of runs of MPCV
#' @param stop.criterion
#' @param max.iter
#' @param max.dim
#' @param method method to be used to determine best run. Possible values are 'h' and 'r'
#' @param scale Should data be scaled?
#' @return a list consisting of
#' \item{segmentation}{of points to clusters}
#' \item{BIC}{Value of \code{\link{myBIC}} criterion}
MPCV.reps <- function(X, numb.Clusters=2, numb.runs=20, stop.criterion=1, max.iter=20, initial.segmentation=NULL,
                      max.dim=1, method=c("likelihood", "singular", "residual"), scale=T, verbose=F){
  if(verbose){
    # create progress bar
    pb <- txtProgressBar(min = 0, max = numb.runs, style = 3)
  }
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
  foreach(i=1:numb.runs) %dopar% {
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
      Hs[i] = H #sum.explained.variance(dane, current.segmentation, max.dim, numb.Clusters)
    }
    
    if(method=="residual"){
      R <- 0 #residual variance
      for (j in 1:numb.Clusters){
        if(is.matrix(dane[,current.segmentation==j]) && dim(dane[,current.segmentation==j])[2]>max.dim){
          a <- summary(prcomp(dane[,current.segmentation==j]))
          R <- R + sum(lm(dane[,current.segmentation==j] ~ current.pcas[[j]])$residuals^2)
        } #otherwise we have a perfect fit
      }
      Res[i] = R #sum.residuals(dane, current.segmentation, max.dim, numb.Clusters)
    }
    
    BICs[i] <- myBIC(dane, current.segmentation, max.dim, numb.Clusters)
    segmentations[[i]] <- current.segmentation
    if(verbose) setTxtProgressBar(pb, i)
  }
  if(verbose) close(pb)
  switch(method,
    singular   = return(list(segmentation = segmentations[[which.max(Hs)]],   BIC = BICs[which.max(Hs)])),
    residual   = return(list(segmentation = segmentations[[which.min(Res)]],  BIC = BICs[which.min(Res)])),
    likelihood = return(list(segmentation = segmentations[[which.max(BICs)]], BIC = BICs[which.max(BICs)])))
}

sum.explained.variance <- function(dane, current.segmentation, max.dim, numb.Clusters){
  H <- 0
  for (j in 1:numb.Clusters){
    if(is.matrix(dane[,current.segmentation==j]) && dim(dane[,current.segmentation==j])[2]>max.dim){
      a <- summary(prcomp(dane[,current.segmentation==j]))
      H <- H + a$importance[3, max.dim]
    }
    else{
      H <- H+1 #all variance in that cluster is explained
    }
  }
  return(H)
}

sum.residuals <- function(dane, current.segmentation, max.dim, numb.Clusters){
  R <- 0
  for (j in 1:numb.Clusters){
    if(is.matrix(dane[,current.segmentation==j]) && dim(dane[,current.segmentation==j])[2]>max.dim){
      a <- summary(prcomp(dane[,current.segmentation==j]))
      R <- R + sum(lm(dane[,current.segmentation==j] ~ current.pcas[[j]])$residuals^2)
    } #otherwise we have a perfect fit
  }
  return(R)
}