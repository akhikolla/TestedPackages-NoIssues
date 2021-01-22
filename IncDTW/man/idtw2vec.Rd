\name{idtw2vec}
\alias{idtw2vec}
\alias{idtw2vec_cm}
\alias{idtw2vec_univ}
\alias{idtw2vec_multiv}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Incremental vector-based DTW 
}
\description{
Update the DTW distance for new observations of two time series. 
}
\usage{
idtw2vec(Q, newObs, dist_method = c("norm1", "norm2", "norm2_square"),
         step_pattern = c("symmetric2", "symmetric1"),
         gcm_lc = NULL, gcm_lr = NULL, nC = NULL, ws = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Q}{Either \code{Q} is (a) a time series (vector or matrix for multivariate time series) or (b) \code{Q} is a cost matrix, so a matrix storing the local distances of the time series \code{Q} and \code{newObs}. If \code{Q} and \code{newObs} are matrices, they need to have the same number of columns. If \code{Q} is a cost matrix, see details...
}
\item{newObs}{time series as vector or matrix, or if \code{Q} is a cost matrix newObs must equals \code{"cm"}. If \code{newObs} is a time series, see details... 
}
\item{dist_method}{character, describes the method of distance measure. See also \code{\link{dtw}}.}
\item{step_pattern}{character, describes the step pattern. See also \code{\link{dtw}}.
}
\item{gcm_lc}{vector, last column of global cost matrix of previous calculation. If NULL (necessary for the initial calculation), then DTW is calculated and the last column and last row are returned to start upcoming incremental calculations. (default = NULL)
}
\item{gcm_lr}{vector, last row of global cost matrix of previous calculation (default = NULL). 
}
\item{nC}{
integer, is the length of the original time series C, of which newObs are the new observations. Length of time series C exclusive new observations, such that \code{length(c(C, newObs)) = nC + length(newObs)}. Necessary if ws is not NULL.  (default = NULL)
}
\item{ws}{
integer, describes the window size for the sakoe chiba window. If NULL, then no window is applied. (default = NULL)
}
}


\details{
If new observations are recorded only for C and the only interest is a fast update of the DTW distance, the last row is not required, neither for the current nor for future incremental calculations.

If \code{Q} is a cost matrix, it needs to store either the distances of \code{Q} and new observations of \code{C} (running calculations, in that case gcm_lc != NULL), or it stores the distances of \code{Q} and the entire time series \code{C} (initial calculation, in that case gcm_lc = NULL).

If \code{newObs} is a time series, it stores either new Observations of \code{C} (running calculations) or the complete time series \code{C} (initial calculation).

no matrices are allocated, no matrices are returned
}
\value{
\item{distance}{the DTW distance}
\item{gcm_lc_new}{the last column of the new global cost matrix}
\item{gcm_lr_new}{the last row of the new global cost matrix. Only if the input vector \code{gcm_lr} is not NUll and represents the last row of the previous global cost matrix, \code{gcm_lr_new} actually is the last row of the updated global cost matrix. Otherwise, if \code{gcm_lr} is NULL then \code{gcm_lr_new} is only the last row of the new part (concerning the new observations) of the global cost matrix.  
}
\item{normalized_distance}{the normalized DTW distance, see also \code{\link{dtw}}}
}
\references{
Sakoe, H.; Chiba, S., Dynamic programming algorithm optimization for spoken word recognition, Acoustics, Speech, and Signal Processing [see also IEEE Transactions on Signal Processing], IEEE Transactions on , vol.26, no.1, pp. 43-49, Feb 1978. http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055
}


\examples{

#--- Do the vector-based incremental DTW 
# calculation and compare it with the basic
Q <- cumsum(rnorm(100))
C <- Q[11:100] + rnorm(90, 0, 0.5)

# initial calculation
res0 <- idtw2vec(Q = Q, newObs = C, gcm_lc = NULL)

# incremental calculation for new observations
nobs <- rnorm(10)
res1 <- idtw2vec(Q, newObs = nobs, gcm_lc = res0$gcm_lc_new)

# compare with result from scratch
res2 <- dtw2vec(Q, c(C, nobs))
res1$distance - res2$distance



#--- Perform an incremental DTW calculation with a 
#  customized distance function. 
d_cos <- function(x, y){
   1 - sum(x * y)/(sqrt(sum(x^2)) * sqrt(sum(y^2))) 
}

x <- matrix(rnorm(100), ncol = 5, nrow = 20)
y <- matrix(rnorm(150), ncol = 5, nrow = 30)
cm1 <- cm(x, y, dist_method = d_cos)

# initial calculation
res0 <- idtw2vec(Q = cm(x, y[1:20,], dist_method = d_cos), 
                 newObs =  "cm")

# incremental calculation for new observations
res1 <- idtw2vec(Q = cm(x, y[21:30,], d_cos), newObs =  "cm", 
         gcm_lc = res0$gcm_lc_new)$distance

# compare with result from scratch
res2 <- dtw2vec(Q = cm1, C = "cm")$distance
res1 - res2

}

\keyword{ classif }% use one of  RShowDoc("KEYWORDS")
\keyword{ cluster }% use one of  RShowDoc("KEYWORDS")
\keyword{ ts }% __ONLY ONE__ keyword per line

