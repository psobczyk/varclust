#' Multiple Latent Components Clustering - kmeans algorithm
#'
#' Performs k-means based subspace clustering. Center of each cluster is some number 
#' of principal components. Similarity measure is calculated using BIC.
#'
#' @param X a matrix with only continuous variables
#' @param number.clusters an integer, number of clusters to be used
#' @param stop.criterion an integer indicating how many changes in partitions triggers stopping the algorithm
#' @param max.iter an integer, maximum number of iterations of k-means
#' @param max.subspace.dim an integer, maximum dimension of subspaces
#' @param initial.segmentation a vector, initial segmentation of variables to clusters
#' @param estimate.dimensions a boolean, if TRUE (value set by default) subspaces dimensions are estimated
#' @param mode a string, possible values : random, kmeans++, sPCA, determines the intialization mode of the algorithm
#' @param show.warnings a boolean - if set to TRUE all warnings are displayed, default value is FALSE
#' @export
#' @return A list consisting of:
#' \item{segmentation}{a vector containing the partition of the variables}
#' \item{pcas}{a list of matrices, basis vectors for each cluster (subspace)}
mlcc.kmeans <- function(X, number.clusters=2, stop.criterion=1, max.iter=40, max.subspace.dim=4, 
                        initial.segmentation=NULL, estimate.dimensions=FALSE, mode = "random", show.warnings = FALSE){
  numbVars = dim(X)[2]
  rowNumb = dim(X)[1]
  pcas <- list(NULL)
  
  if(is.null(initial.segmentation)){
    switch(mode,
           "random"={
             los = sample(1:numbVars,number.clusters)
             pcas = lapply(1:number.clusters, function(i) matrix(X[,los[i]], nrow=rowNumb))
           },
           "kmeans++"={
               los = sample(1:numbVars,1)
               pcas = list(matrix(X[,los], nrow = rowNumb))
               for( k in 2:number.clusters ){
                 dists = vapply(1:numbVars, function(i) calculate.distance.kmeanspp(X[,i],pcas,k-1), 0.9)
                 distsum = sum(dists)
                 new_center_index = sample(1:numbVars, 1 , prob = dists/distsum)
                 pcas[[k]] <- matrix(X[,new_center_index], nrow = rowNumb)
               }
             },
             "sPCA"={
                if(rowNumb < number.clusters){
                  stop("Number of cluster cannot be bigger than number of rows")
                }
                count <- min(round(number.clusters*log(numbVars)), rowNumb)
                PCs <- SPC(x=X, sumabsv=sqrt(numbVars), niter = 50, K = count, trace = FALSE)$u
                pcas = list(matrix(PCs[,1],nrow = rowNumb))
                for( k in 2:number.clusters ){
                  dists = vapply(1:numbVars, function(i) calculate.distance.kmeanspp(X[,i],pcas,k-1), 0.9)
                  new_center_index = which.max(dists)
                  pcas[[k]] <- matrix(X[,new_center_index], nrow = rowNumb)
                }
            })
    segmentation <- sapply(1:numbVars,  function(j) choose.cluster.BIC(X[,j],pcas, number.clusters, show.warnings))
  }
  else{
    segmentation=initial.segmentation
  }
  
  new.segmentation <- segmentation
  for (iter in 1:max.iter){
    pcas <- lapply(1:number.clusters, function(k){
      Xk = X[,segmentation==k, drop=F]
      sub.dim <- dim(Xk)
      if(sub.dim[2] > 0){
        a <- summary(prcomp(x=Xk))
        if (estimate.dimensions) {
            max.dim <- min(max.subspace.dim, floor(sqrt(sub.dim[2])), sub.dim[1])
            cut <- pesel(X = Xk, npc.min = 1, npc.max = max.dim, scale = FALSE, method = "heterogenous")$nPCs
        }
        else {
          cut <- min(max.subspace.dim, floor(sqrt(sub.dim[2])), sub.dim[1])
        }
        return(matrix(a$x[,1:cut], nrow=rowNumb))
      }
      else{
        return(matrix(rnorm(rowNumb*dim), nrow = rowNumb, ncol = 1))
      }
    })
    new.segmentation <- sapply(1:numbVars, function(j) choose.cluster.BIC(X[,j], pcas, number.clusters, show.warnings))
    if(sum(new.segmentation!=segmentation)<stop.criterion) break
    segmentation = new.segmentation
  }
  return(list(segmentation=segmentation, 
              pcas=pcas))
}