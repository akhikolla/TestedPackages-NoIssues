\name{dtw_dismat}
\alias{dtw_dismat}
\alias{dtw_disvec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DTW Distance Matrix/ Distance Vector
}
\description{
Calculate a matrix of pairwise DTW distances for a set of univariate or multivariate time series. The output matrix (or dist object) of DTW distances can easily be applied for clustering the set of time series. Or calculate a vector of DTW distances of a set of time series all relative to one query time series. Parallel computations are possible.}
\usage{
dtw_dismat(lot, dist_method = c("norm1", "norm2", "norm2_square"),
         step_pattern = c("symmetric2", "symmetric1"), normalize = TRUE,
         ws = NULL, threshold = NULL,
         return_matrix = FALSE, ncores = NULL, useRcppParallel = TRUE)

dtw_disvec(Q, lot, dist_method = c("norm1", "norm2", "norm2_square"),
         step_pattern = c("symmetric2", "symmetric1"), normalize = TRUE,
         ws = NULL, threshold = NULL, ncores = NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Q}{time series, vector (univariate) or matrix (multivariate)}
\item{lot}{List of time series. Each entry of the list is a time series as described in \code{\link{dtw2vec}}.}
\item{dist_method}{character, describes the method of distance measure. See also \code{\link{dtw}}.}
\item{step_pattern}{character, describes the step pattern. See also \code{\link{dtw}}.}
\item{normalize}{logical, whether to return normalized pairwise distances or not. If \code{step_pattern == 'symmetric1'} only non-normalized distances can be returned (default = TRUE)}
\item{ws}{
integer, describes the window size for the sakoe chiba window. If NULL, then no window is applied. (default = NULL)
}
\item{threshold}{numeric, the threshold for early abandoning. In the calculation of the global cost matrix a possible path stops as soon as the threshold is reached. Facilitates faster calculations in case of low threshold. (default = FALSE)
}
\item{return_matrix}{logical, If FALSE (default) the distance matrix is returned as \code{dist} object. If TRUE a symmetric matrix of differences is returned.  
}
\item{ncores}{integer, number of cores to be used for parallel computation of the distance matrix. If \code{ncores = NULL} (default) then ncores is set to the number of available cores minus 1. If \code{ncores = 0} then no parallel computation is performed and standard \code{sapply} instead of \code{parallel::parSapply} is applied.}
\item{useRcppParallel}{logical, if the package RcppParallel (TRUE, default) or parallel  (FALSE) is used for parallel computation}
}


\details{
By setting the parameter \code{return_matrix = FALSE} (default) the output value \code{dismat} of \code{dtw_dismat} is a \code{dist} object and can easily be passed to standard clustering functions (see examples).

No matrices are allocated for calculating the pairwise distances.
}
\value{
\item{input}{the function input parameters}
\item{dismat}{the matrix of pairwise DTW distances, either as matrix or dist object}
\item{disvec}{the vector DTW distances}
}



\examples{
\dontrun{

#--- Example for clustering a set of time series by feeding well known 
# clustering methods with DTW-distance objects. First we simulate 
# two prototype random walks and some cluster members. The cluster
# members are simulated by adding noise and randomly stretching and 
# comressing the time series, to get time warped time series of 
# varying lengths. The built clusters are 1:6 and 7:12.
set.seed(123)
N <- 100
rw_a <- cumsum(rnorm(N))
rw_b <- cumsum(rnorm(N))
sth <- sample(seq(0, 0.2, 0.01), size = 10)
cmp <- sample(seq(0, 0.2, 0.01), size = 10)
lot <- c(list(rw_a), 
         lapply(1:5, function(i){ 
           simulate_timewarp(rw_a + rnorm(N), sth[i], cmp[i]) 
         }),
         list(rw_b),
         lapply(6:10, function(i){ 
           simulate_timewarp(rw_b + rnorm(N), sth[i], cmp[i]) 
         }))


# Next get the distance matrix, as dist object. Per default all 
# minus 1 available cores are used:
result <- dtw_dismat(lot = lot, dist_method = "norm2", ws = 50, 
                     return_matrix = FALSE)
class(result$dismat)


# Finally you can cluster the result with the following
#  well known methods:
require(cluster)
myclus <- hclust(result$dismat)
plot(myclus)
summary(myclus)

myclus <- agnes(result$dismat)
plot(myclus)
summary(myclus)

myclus <- pam(result$dismat, k=2)
plot(myclus)
summary(myclus)
myclus$medoids


}
}

\keyword{ classif }% use one of  RShowDoc("KEYWORDS")
\keyword{ cluster }% use one of  RShowDoc("KEYWORDS")
\keyword{ ts }% __ONLY ONE__ keyword per line

