## ----, results='hide'----------------------------------------------------
library(VarClust)
library(mclust)

## ------------------------------------------------------------------------
comp <- read.table("http://factominer.free.fr/docs/gene.csv",sep=";",header=T,row.names=1) 
benchamrkClustering <- c(rep(1, 68), rep(2, 356))    
comp <- data.frame(comp[,-ncol(comp)])    
mlcc.fit <- MPCV.BIC(comp, numb.Clusters = 1:3, numb.runs = 100, max.dim = 2, numbCores=4, greedy=F, estimateDimensions=F)


