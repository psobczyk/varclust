\name{MPCV}
\alias{MPCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs subspace clustering (variables clustering)
}
\description{
clusters variables based on k-means algorithm where distance is measured with R^2 (from least squares)
}
\usage{
MPCV(X, numberClusters = 2, stop_criterion = 2, var_threshold = 0.5, max_iter = 20, maxSubspaceDim = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
data
}
  \item{numberClusters}{
number of clusters to be used (needed as functions uses k-means algorithm)
}
  \item{stop_criterion}{
if minimal change in segmentation reaches stop_criterion k-means stops 
}
  \item{var_threshold}{
threshold used to determine the dimension of subspace. Percent of variance that needs to be explained.
}
  \item{max_iter}{
maximum number of iterations for k-means
}
  \item{maxSubspaceDim}{
maximum dimension of subspaces
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{segmentation }{Segmentation of variables into clusters}
\item{pcas }{PCA components of different clusters}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Piotr Sobczyk
}
\note{
Be cautious...
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library(MASS)

n = 50;
sigma = 0.05;

K = 2; #Number of subspaces
dims = rep(1, K)
N = 100/K;

X = NULL;
s = NULL
for (j in 1:K){
  Z = matrix(mvrnorm(dims[j], mu=rep(0,n), Sigma=diag(rep(1,n))), ncol=dims[j])
  X = cbind(X, Z %*% mvrnorm(dims[j],rep(0,N),diag(rep(1,N))) + mvrnorm(n,rep(0,N),diag(rep(sigma,N))))
  s = c(s, rep(j, N));
}

group <- MPCV(X, numberClusters=K)[[1]]
missclassify_heuristic(group, N, K)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ cluster }
