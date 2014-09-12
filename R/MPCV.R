#' Subspace clustering based on multiple principal components
#'
#' Performs subspace clustering
#'
#' @param X data
#' @param numberClusters number of clusters to be used
#' @param stop.criterion
#' @param max.iter
#' @param maxSubspaceDim
#' @param initial.segmentation
#' @param estimateDimensions
#' @param known.sigma
#' @return a list consisting of
#' \item{segmentation}{of points to clusters}
#' \item{pcas}{basis of each cluster}
MPCV <- function(X, numberClusters = 2, stop.criterion  = 1, max.iter = 40, maxSubspaceDim=4, 
                 initial.segmentation=NULL, estimateDimensions=FALSE, known.sigma=1){
  X = scale(X)
  numbVars = dim(X)[2]
  rowNumb = dim(X)[1]
  pcas <- list(NULL)
  
  if(is.null(initial.segmentation)){ #if there is no initial segmentation
    los = sample(1:numbVars,numberClusters)
    pcas = lapply(1:numberClusters, function(i) matrix(X[,los[i]], nrow=rowNumb))
    #all subspaces are univariate so residuals are as good as BIC
    segmentation <- sapply(1:numbVars,  function(j) choose.cluster(X[,j],pcas, numberClusters)) 
  }
  else
    segmentation=initial.segmentation
  
  new.segmentation <- segmentation #auxillary variable
  
  for (iter in 1:max.iter){
    for (i in 1:numberClusters){
      if(is.matrix(X[,segmentation==i]) && (dim(X[,segmentation==i])[2]>=maxSubspaceDim)){
        a <- summary(prcomp(x=X[,segmentation==i]))
        if(estimateDimensions){
          cut <- which.max(lapply(1:maxSubspaceDim, 
                                  function(x) myBIC(X[,segmentation==i], rep(1,sum(segmentation==i)), max.dim=x, numb.clusters=1, sigma=known.sigma)))
        }
        else{
          cut <- maxSubspaceDim
        }
        pcas[[i]] = matrix(a$x[,1:cut], nrow=rowNumb)
      }
      else{
        m <- as.numeric(names(sort(table(new.segmentation),decreasing=T)))[1] #most populous cluster will be cut into two
        try({
          podzial <- kmeansvar(X[,new.segmentation==m],init=2)$cluster
          new.segmentation[which(new.segmentation==m)[podzial==2]] = i
          segmentation[which(new.segmentation==m)[podzial==2]] = i
        })
        if(length(X[,segmentation==i])<maxSubspaceDim*rowNumb) #if that still doesn't help new cluster is choosen at random
          pcas[[i]] = matrix(rnorm(rowNumb*maxSubspaceDim), nrow=rowNumb)
        else{
          a <- summary(prcomp(x=X[,segmentation==i]))
          if(estimateDimensions){
            cut <- which.max(lapply(1:maxSubspaceDim, 
                                    function(x) myBIC(X[,segmentation==i], rep(1,sum(segmentation==i)), max.dim=x,numb.clusters=1), sigma=known.sigma))
          }
          else{
            cut <- maxSubspaceDim
          }
          pcas[[i]] = matrix(a$x[,1:cut], nrow=rowNumb)
        } 
      }
    }
    if(estimateDimensions){ #if we try to estimate dimension we use BIC
      new.segmentation <- sapply(1:numbVars, function(j) choose.cluster.BIC(X[,j], pcas, numberClusters, known.sigma))
    }
    else{ #otherwise we use R^2
      new.segmentation <- sapply(1:numbVars, function(j) choose.cluster(X[,j], pcas, numberClusters))
    }
    if(sum(new.segmentation!=segmentation)<stop.criterion) break
    segmentation = new.segmentation
  }
  return(list(segmentation=segmentation, 
              pcas=pcas))
}