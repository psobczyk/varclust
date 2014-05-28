
MPCV.reps <- function(X, true.segmentation, numb.Clusters=2, numb.runs=20, stop.criterion=1, max.iter=20, max.dim=1, method='h', scale=T){
  if(scale){
    dane = scale(X)
  }
  else dane= X
  #goodness of clustering 
  missclassifications <- NULL
  adjs <- NULL
  #goodness of fit
  Hs <- NULL
  Res <- NULL
  BICs <- NULL 
  
  segmentations <- NULL
  for(i in 1:numb.runs){
    H <- 0 #explained variance
    R <- 0 #residual variance
    
    MPCV.res <- MPCV(X=dane, numberClusters=numb.Clusters, maxSubspaceDim=max.dim)
    current.segmentation <- MPCV.res$segmentation
    current.pcas <- MPCV.res$pcas
   
    for (j in 1:numb.Clusters){
      if(is.matrix(dane[,current.segmentation==j]) && dim(dane[,current.segmentation==j])[2]>max.dim){
        a <- summary(prcomp(dane[,current.segmentation==j]))
        H <- H + a$importance[3,max.dim]
      }
      else{
        H <- H+1 #all variance in that cluster is explained
      }
    }
    
    for (j in 1:numb.Clusters){
      if(is.matrix(dane[,current.segmentation==j]) && dim(dane[,current.segmentation==j])[2]>max.dim){
        a <- summary(prcomp(dane[,current.segmentation==j]))
        R <- R + sum(lm(dane[,current.segmentation==j] ~ current.pcas[[j]])$residuals^2)
      } #otherwise we have a perfect fit
    }
    
    Hs[i] = H
    Res[i] = R 
    missclassifications[i] <- missclassify.heuristic(current.segmentation, dim(dane)[2]/numb.Clusters, n=numb.Clusters)
    adjs[i] <- adjustedRandIndex(current.segmentation, true.segmentation)
    BICs[i] <- myBIC(dane, current.segmentation,max.dim, numb.Clusters)
    segmentations[[i]] <- current.segmentation
    if(i%%10==0) print(paste("Done ", i))
  }  
  if(method=='h'){
    return(list(segmentation = segmentations[[which.max(Hs)]],
                misclassification = missclassifications[which.max(Hs)],
                ADI = adjs[which.max(Hs)],
                BIC = BICs[which.max(Hs)]))
  }
  else{
    return(list(segmentation = segmentations[[which.min(Res)]],
                misclassification = missclassifications[which.min(Res)],
                ADI = adjs[which.min(Res)],
                BIC = BICs[which.min(Res)]))
  }
}