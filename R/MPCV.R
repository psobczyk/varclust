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
#' @return a list consisting of
#' \item{segmentation}{of points to clusters}
#' \item{pcas}{basis of each cluster}
MPCV <- function(X, numberClusters = 2, stop.criterion  = 1, max.iter = 40, maxSubspaceDim=4, initial.segmentation=NULL){
  X = scale(X)
  numbVars = dim(X)[2]
  rowNumb = dim(X)[1]
  pcas <- list(NULL)
  
  if(is.null(initial.segmentation)){
    los = sample(1:numbVars,numberClusters)
    pcas = lapply(1:numberClusters, function(i) matrix(X[,los[i]], nrow=rowNumb))
    segmentation <- sapply(1:numbVars,  function(j) choose.cluster(X[,j],pcas, numberClusters))
  }
  else
    segmentation=initial.segmentation
  
  new.segmentation <- segmentation
  
  for (iter in 1:max.iter){
    pcas <- lapply(1:numberClusters, function(i){
      if(is.matrix(X[,segmentation==i]) && (dim(X[,segmentation==i])[2]>=maxSubspaceDim)){
        a <- summary(prcomp(x=X[,segmentation==i]))
        cut <- maxSubspaceDim
        return(matrix(a$x[,1:cut], nrow=rowNumb))
      }
      else{
        m <- as.numeric(names(sort(table(new.segmentation),decreasing=T)))[1] #most populous cluster
        try({
          podzial <- kmeansvar(X[,new.segmentation==m],init=2)$cluster
          new.segmentation[which(new.segmentation==m)[podzial==2]] = i
          segmentation[which(new.segmentation==m)[podzial==2]] = i
        })
        if(length(X[,segmentation==i])<maxSubspaceDim*rowNumb)
          return(matrix(rnorm(rowNumb*maxSubspaceDim), nrow=rowNumb))
        else{
          a <- summary(prcomp(x=X[,segmentation==i]))
          cut <- maxSubspaceDim
          return(matrix(a$x[,1:cut], nrow=rowNumb))
        } 
      }
    })
    new.segmentation <- sapply(1:numbVars, function(j) choose.cluster(X[,j],pcas, numberClusters))
    if(sum(new.segmentation!=segmentation)<stop.criterion) break
    segmentation = new.segmentation
  }
  return(list(segmentation=segmentation, 
              pcas=pcas))
}