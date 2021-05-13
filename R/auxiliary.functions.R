#' mBIC for subspace clustering
#'
#' Computes the value of modified Bayesian Information Criterion (mBIC) for
#' given data set partition and clusters' dimensionalities. In each cluster we
#' assume that variables are spanned by few factors. Considering maximum
#' likelihood we get that those factors are in fact principal components.
#' Additionally, it uses by default an informative prior distribution on models.
#'
#'
#' @param X A matrix with only quantitative variables.
#' @param segmentation A vector, segmentation for which likelihood is computed.
#'   Clusters numbers should be from range [1, numb.clusters].
#' @param dims A vector of integers, dimensions of subspaces. Number of
#'   principal components (fixed or chosen by PESEL criterion) that span each
#'   subspace.
#' @param numb.clusters An integer, number of clusters.
#' @param max.dim An integer, upper bound for allowed dimension of a subspace.
#' @param flat.prior A boolean, if TRUE (default is FALSE) then flat prior on
#'   models is used.
#' @keywords internal
#' @return Value of mBIC
cluster.pca.BIC <- function(X, segmentation, dims, numb.clusters, max.dim, flat.prior = FALSE) {
  if (!is.matrix(X)) {
    # if X is one variable it is stored as vector
    X <- matrix(X, ncol = 1)
  }
  D <- dim(X)[1]
  p <- dim(X)[2]

  formula <- rep(0, numb.clusters)
  for (k in 1:numb.clusters) {
    # one cluster
    dimk <- dims[k]
    Xk <- X[, segmentation == k, drop = F]
    if (dim(Xk)[2] > dimk && dim(Xk)[1] > dimk) {
      formula[k] <- pesel(
        X = Xk, npc.min = dimk, npc.max = dimk, scale = FALSE,
        method = "heterogenous"
      )$vals[1]
    } else {
      warning("The dimensionality of the cluster was greater or equal than
              max(number of observation, number of variables) in the cluster.
              Ignoring the cluster during mBIC calculation")
      formula[k] <- 0
    }
  }
  # apriori
  apriori.segmentations <- -p * log(numb.clusters)
  apriori.dimensions <- -log(max.dim) * numb.clusters
  BIC <- sum(formula)
  if (!flat.prior) {
    BIC <- BIC + apriori.segmentations + apriori.dimensions
  }
  return(BIC)
}

#' Chooses a subspace for a variable
#'
#' Selects a subspace closest to a given variable. To select the subspace, the method
#' considers (for every subspace) a subset of its principal components and tries
#' to fit a linear model with the variable as the response. Then the method chooses
#' the subspace for which the value of BIC was the highest.
#'
#' @param variable A variable to be assigned.
#' @param pcas Orthogonal basis for each of the subspaces.
#' @param number.clusters Number of subspaces (clusters).
#' @param show.warnings A boolean - if set to TRUE all warnings are displayed, default value is FALSE.
#' @param common_sigma A boolean - if set to FALSE, seperate sigma is estimated for each cluster,
#' default value is TRUE
#' @keywords internal
#' @return index Number of most similar subspace to variable.
choose.cluster.BIC <- function(variable, pcas, number.clusters, show.warnings = FALSE, common_sigma = TRUE) {
  BICs <- NULL
  if (common_sigma) {
    res <- fastLmPure(cbind(1, as.matrix(Matrix::bdiag(pcas))), rep(variable, number.clusters), method = 0L)$residuals
    n <- length(variable)
    sigma.hat <- sqrt(sum(res^2) / (n * number.clusters))
    if (sigma.hat < 1e-15 && show.warnings) {
      warning("In function choose.cluster.BIC: estimated value of noise in cluster is <1e-15. It might corrupt the result.")
    }
    for (i in 1:number.clusters) {
      nparams <- ncol(pcas[[i]])
      res.part <- res[((i - 1) * n + 1):(i * n)]
      loglik <- sum(dnorm(res.part, 0, sigma.hat, log = T))
      BICs[i] <- loglik - nparams * log(n) / 2
    }
  } else {
    for (i in 1:number.clusters) {
      nparams <- ncol(pcas[[i]])
      n <- length(variable)
      res <- fastLmPure(pcas[[i]], variable, method = 0L)$residuals
      sigma.hat <- sqrt(sum(res^2) / n)
      if (sigma.hat < 1e-15 && show.warnings) {
        warning("In function choose.cluster.BIC: estimated value of noise in cluster is <1e-15. It might corrupt the result.")
      }
      loglik <- sum(dnorm(res, 0, sigma.hat, log = T))
      BICs[i] <- loglik - nparams * log(n) / 2
    }
  }
  return(which.max(BICs))
}

#' Calculates principal components for every cluster
#'
#' For given segmentation this function estimates dimensionality of each cluster (or chooses fixed dimensionality)
#' and for each cluster calculates the number of principal components equal to the this dimensionality
#'
#' @param X A data matrix.
#' @param segmentation A vector, segmentation of variables into clusters.
#' @param number.clusters An integer, number of subspaces (clusters).
#' @param max.subspace.dim An integer, upper bound for allowed dimension of subspace.
#' @param estimate.dimensions A boolean, if TRUE subspaces dimensions are estimated using PESEL.
#' @keywords internal
#' @return A subset of principal components for every cluster.
calculate.pcas <- function(X, segmentation, number.clusters, max.subspace.dim, estimate.dimensions) {
  rowNumb <- dim(X)[1]
  pcas <- lapply(1:number.clusters, function(k) {
    Xk <- X[, segmentation == k, drop = F]
    sub.dim <- dim(Xk)
    if (sub.dim[2] > 0) {
      a <- summary(prcomp(x = Xk))
      if (estimate.dimensions) {
        max.dim <- min(max.subspace.dim, floor(sqrt(sub.dim[2])), sub.dim[1])
        cut <- max(0, pesel(
          X = Xk, npc.min = 0, npc.max = max.dim, scale = FALSE,
          method = "heterogenous"
        )$nPCs)
      } else {
        cut <- min(max.subspace.dim, floor(sqrt(sub.dim[2])), sub.dim[1])
      }
      return(matrix(a$x[, seq_len(cut)], nrow = rowNumb))
    } else { #if there are no variables initiate a cluster at random
      return(matrix(rnorm(rowNumb), nrow = rowNumb, ncol = 1))
    }
  })
  return(pcas)
}

#' Plot mlcc.fit class object
#'
#' @param x mlcc.fit class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
plot.mlcc.fit <- function(x, ...) {
  clusterNumbs <- lapply(x$all.fit, function(y) y$nClusters)
  BICVals <- lapply(x$all.fit, function(y) y$BIC)
  plot.default(clusterNumbs, BICVals, type = "b", xaxt = "n", ylab = "BIC", xlab = "Number of clusters")
  axis(side = 1, labels = clusterNumbs, at = clusterNumbs)
}

#' Print mlcc.fit class object
#'
#' @param x mlcc.fit class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
print.mlcc.fit <- function(x, ...) {
  cat("$nClusters: ", x$nClusters, "\n")
  cat("$segmentation:\n")
  print(x$segmentation)
  cat("$BIC: ", x$BIC, "\n")
  cat("$subspacesDimensions:\n", unlist(x$subspacesDimensions), "\n")
}

#' Print mlcc.reps.fit class object
#'
#' @param x mlcc.reps.fit class object
#' @param ... Further arguments to be passed to or from other methods. They are
#'   ignored in this function.
#' @export
#' @keywords internal
print.mlcc.reps.fit <- function(x, ...) {
  cat("$segmentation:\n")
  print(x$segmentation)
  cat("$BIC: ", x$BIC, "\n")
  cat("$basis:\n")
  cat(str(x$basis))
}

#' Print clusters obtained from MLCC
#'
#' @param data The original data set.
#' @param segmentation A vector, segmentation of variables into clusters.
#' @export
show.clusters <- function(data, segmentation) {
  data <- as.data.frame(data)
  max_cluster_size <- max(as.data.frame(table(segmentation))$Freq)
  clusters <- lapply(1:max(segmentation), function(i) {
    colnames_in_cluster <- colnames(data)[segmentation == i]
    current_cluster_size <- length(colnames_in_cluster)
    c(
      colnames_in_cluster,
      rep("-", times = max_cluster_size - current_cluster_size)
    )
  })
  clusters <- as.data.frame(clusters)
  colnames(clusters) <- paste("cluster", 1:max(segmentation), sep = "_")
  print(clusters)
}
