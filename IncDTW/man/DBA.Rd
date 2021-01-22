\name{dba}
\alias{dba}
\alias{DBA}
\alias{centroid}
\alias{print.dba}
\alias{summary.dba}
\alias{is.dba}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dynamic Time Warping Barycenter Averaging
}
\description{
Average multiple time series that are non-linearly aligned by Dynamic Time Warping. Find the centroid of a list of time series.
}
\usage{
dba(lot, m0 = NULL, iterMax = 10, eps = NULL, 
                dist_method = c("norm1", "norm2", "norm2_square"),
                step_pattern = c("symmetric2", "symmetric1"),
                ws = NULL,
                iter_dist_method = c("dtw_norm1", "dtw_norm2",
                                     "norm1","norm2", "max", "min"), 
                plotit = FALSE)
                
# deprecated
DBA(lot, m0 = NULL, iterMax = 10, eps = NULL, 
                dist_method = c("norm1", "norm2", "norm2_square"),
                step_pattern = c("symmetric2", "symmetric1"),
                ws = NULL,
                iter_dist_method = c("dtw_norm1", "dtw_norm2",
                                     "norm1","norm2", "max", "min"), 
                plotit = FALSE)
                
centroid(lot, dist_method = c("norm1", "norm2", "norm2_square"),
              step_pattern = c("symmetric2", "symmetric1"),
              normalize = TRUE, ws = NULL, ncores = NULL,
              useRcppParallel = TRUE)


              
\method{print}{dba}(x, digits = getOption("digits"),  ...)

\method{summary}{dba}(object, ...)              

is.dba(x)             
              
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{lot}{List of time series. Each entry of the list is a time series as described in \code{\link{dtw2vec}}.}
\item{m0}{time series as vector or matrix. If \code{m0} is \code{NULL}, the initial time series \code{m0} is determined by \code{\link{centroid}} as the centroid of \code{lot}, which is the one time series of \code{lot} with the minimum average DTW distance to all other time series of \code{lot}.}
\item{iterMax}{integer, number of maximum iterations}
\item{eps}{numeric, threshold parameter that causes the algorithm to break if the distance of two consecutive barycenters are closer than eps}
\item{dist_method}{character, describes the method of distance measure. See also \code{\link{dtw}}.}
\item{step_pattern}{character, describes the step pattern. See also \code{\link{dtw}}.}
\item{ws}{
integer, describes the window size for the sakoe chiba window. If NULL, then no window is applied. (default = NULL)
}
\item{iter_dist_method}{character, that describes how the distance between two consecutive barycenter iterations are defined (default = "dtw")}
\item{plotit}{logical, if the iterations should be plotted or not (only possible for univariate time series)}
\item{normalize}{logical, default is TRUE, passed to \code{\link{dtw_dismat}}}
\item{ncores}{integer, default = NULL,  passed to \code{\link{dtw_dismat}}}
\item{useRcppParallel}{logical, default is TRUE, passed to \code{\link{dtw_dismat}}}
\item{x}{
output from \code{dba}
}
\item{object}{
any R object
}
\item{digits}{
passed to \code{round} and \code{print}
}
\item{...}{
additional arguments, e.g. passed to \code{print} or \code{summary}
}

}


                   
                   
\details{
The parameter \code{iter_dist_method} describes the method to measure the progress between two iterations. For two consecutive centroid candidates \code{m1} and \code{m2} the following methods are implemented:

'dtw_norm1':
\code{dtw2vec(m1, m2, dist_method = "norm1", step_pattern = "symmetric2")$normalized_distance}


'idm_dtw2':
   \code{dtw2vec(m1, m2, dist_method = "norm2", step_pattern = "symmetric2")$normalized_distance}


'idm_norm1':
   \code{sum(abs(m1-m2))/(ncol(m1) * 2 * nrow(m1))}


'idm_norm2':
   \code{sqrt(sum((m1-m2)^2))/(ncol(m1) * 2 * nrow(m1))}


'idm_max':
   \code{max(abs(m1-m2))}


'idm_min':
   \code{min(abs(m1-m2))}
   
}

\value{
\item{call}{function call}
\item{m1}{new centroid/ bary center of the list of time series}
\item{iterations}{list of time series that are the best centroid of the respective iteration}
\item{iterDist_m2lot}{list of distances of the iterations to lot}
\item{iterDist_m2lot_norm}{list of normalized distances of the iterations to lot}
\item{iterDist_m2m}{vector of distances of the iterations to their ancestors}
\item{centroid_index}{integer giving the index of the centroid time series of \code{lot}}
\item{dismat_result}{list of results of \code{\link{dtw_dismat}} called by \code{centroid}}
}

\references{
Sakoe, H.; Chiba, S., Dynamic programming algorithm optimization for spoken word recognition, Acoustics, Speech, and Signal Processing [see also IEEE Transactions on Signal Processing], IEEE Transactions on , vol.26, no.1, pp. 43-49, Feb 1978. \url{https://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055}

Petitjean, F; Ketterlin, A; Gancarski, P, A global averaging method for dynamic time warping, with applications to clustering, Pattern Recognition, Volume 44, Issue 3, 2011, Pages 678-693, ISSN 0031-3203
}


\examples{
\dontrun{
data("drink_glass")
# initialize with any time series 
m1 <- dba(lot = drink_glass[1:10], m0 = drink_glass[[1]], 
          dist_method = "norm2", iterMax = 20)

# initialize with the centroid

tmp <- centroid(drink_glass)
cent <- drink_glass[[tmp$centroid_index]]
m1 <- dba(lot = drink_glass[1:10], m0 = cent, 
          dist_method = "norm2", iterMax = 20)

# plot all dimensions of the barycenters m_n per iteration:
plot(m1)

# plot the distances of the barycenter of one iteration m_n
# to the barycenter of the previous iteration m_n-1:
plot(m1, type = "m2m")

# plot the average distances of the barycenter m_n 
# to the list of time series:
plot(m1, type = "m2lot")

}
}

\keyword{ classif }% use one of  RShowDoc("KEYWORDS")
\keyword{ cluster }% use one of  RShowDoc("KEYWORDS")
\keyword{ ts }% __ONLY ONE__ keyword per line

