\name{initialize_plane}
\alias{initialize_plane}
\alias{increment}
\alias{increment.planedtw}
\alias{decrement}
\alias{decrement.planedtw}
\alias{refresh}
\alias{refresh.planedtw}
\alias{reverse}
\alias{reverse.planedtw}
\alias{print.planedtw}
\alias{is.planedtw}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Initialize and navigate in the plane of possible fits
}
\description{
Initialize and navigate in the plane of possible fits to detect subsequences (of different lengths) in a long time series that are similar (in terms of DTW distance) to a query pattern: Initialize the plane of possible fits as \code{.planedtw} object. Increment and decrement the time series observations and respective DTW calculation. Reverse the time order to increment or decrement observations at the other end of the time horizon. Refresh the DTW calculation without changing the time series. 
}
\usage{
initialize_plane(Q, C, dist_method = c("norm1", "norm2", "norm2_square"),
                 step_pattern = c("symmetric2", "symmetric1"), ws = NULL)


\method{increment}{planedtw}(x, newObs, direction = c("C", "Q"), ...)

\method{decrement}{planedtw}(x, direction = c("C", "Q", "both"), 
          refresh_dtw = FALSE, nC = NULL, nQ = NULL, ...)

\method{refresh}{planedtw}(x, ...)

\method{reverse}{planedtw}(x, ...)

is.planedtw(x)  

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Q}{
a time series (vector or matrix for multivariate time series) 
}
\item{C}{
a time series (vector or matrix for multivariate time series) 
}

\item{dist_method}{
character, describes the method of distance measure. See also \code{\link{dtw}}.
}
\item{step_pattern}{
character, describes the step pattern. See also \code{\link{dtw}}.
}
\item{ws}{
integer, describes the window size for the sakoe chiba window. If NULL, then no window is applied. (default = NULL)
}
\item{x}{
object of class planedtw (output from \code{initialize_plane})
}
\item{newObs}{
a time series (vector or matrix for multivariate time series). If \code{Q} and \code{C} are vectors, newObs must be a vector. If \code{Q} and \code{C} are matrices with \code{nc} columns, then \code{newObs} must also have \code{nc} columns. See details for the correct time order of \code{newObs}.
}
\item{direction}{
character, gives the direction of increment or decrement. \code{decrement()} is a wrapper for \code{\link{dtw_partial}} and the \code{direction} parameter is translated to the respective \code{partial_Q} and \code{partial_C} parameters.
}
\item{refresh_dtw}{
logical (default = FALSE), after decrementing the time series, should the DTW calculation be refreshed, or not.
}
\item{nC}{
integer, default = NULL, if not NULL, then \code{decrement} subsets the time series \code{C} to the range of \code{1:nC}, drops invalid interim calculation results, and refreshes if \code{refresh_dtw = TRUE}.
}
\item{nQ}{
analog to \code{nC}
}
\item{...}{
additional arguments (currently not used)
}
}


\details{
All functions are wrapper functions for \code{\link{idtw2vec}} and \code{\link{dtw_partial}}. 

\itemize{

\item{
\code{initialize_plane} calculates the DTW distance between Q and C and saves the last column and row of the global cost matrix. It returns an object of class \code{planedtw} that contains all necessary information to incrementally update the DTW calculation with new observations. Also for decrementing the calcultions for skipping some observations at the end.
}
\item{
\code{increment} updates the DTW calculation by appending new observations to C or Q (depends on the parameter \code{direction}) and calculating DTW by recycling previous results represented by \code{gcm_lc_new} and \code{gcm_lr_new}. A wrapper for \code{\link{idtw2vec}}
}
\item{
\code{decrement} is a wrapper for \code{\link{dtw_partial}} and also returns a \code{planedtw} object.
}
\item{
\code{refresh} serves to recalculate the \code{gcm_lc_new} and \code{gcm_lr_new} from scratch, if these objects are \code{NULL} (e.g. after decrementing with \code{refresh_dtw = FALSE}).
}
\item{
\code{reverse} reverses the order of Q and C, and refreshes the calculation for the new order. This is useful for appending observations to Q or C at the other end, the beginning. For incrementing in the reverse order also apply the function \code{increment}. Then the time series in the parameter \code{newObs} also needs to be in reverse order. 
Assent et al. (2009) proved that the DTW distance is reversible for the step pattern "symmetric1", so \code{dtw(Q, C) = dtw(rev(Q), rev(C))}. Also see examples. For the step pattern "symmetric2" DTW is not exactly reversible, but empirical studies showed that the difference is realtive small. For further details please see the appendix A of the vignette "IncDTW: An R Package for Incremental Calculation of Dynamic Time Warping" on  \href{https://CRAN.R-project.org/package=IncDTW}{CRAN}.
}
}


}
\value{
\item{distance}{the DTW distance}
\item{normalized_distance}{the DTW distance devided by the sum of the lengths of Q and C (see also \code{\link{dtw}}).}
\item{gcm_lc_new}{the last column of the new global cost matrix}
\item{gcm_lr_new}{the last row of the new global cost matrix}
\item{Q}{the time series}
\item{C}{the time series}
\item{control}{list of input parameters and the lengths of the time series}
}

\references{
Assent, Ira, et al. "Anticipatory DTW for efficient similarity search in time series databases." Proceedings of the VLDB Endowment 2.1 (2009): 826-837.
}



\examples{
\dontrun{

#--- 1. example: Increment too far and take a step back:
rw <- function(nn) cumsum(rnorm(nn))
Q <- sin(1:100)
C <- Q[1:90] + rnorm(90, 0, 0.1)
WS <- 40


# start with the initial calculation
x <- initialize_plane(Q, C, ws = WS)

# Then the incremental calculation for new observations
y1 <- Q[91:95] + rnorm(5, 0, 0.1)# new observations
x <- increment(x, newObs = y1)

# Again new observations -> just increment x
y2 <- c(Q[96:100] + rnorm(5, 0, 0.1), rw(10))# new observations
x <- increment(x, newObs = y2)

# Compare the results with the calculation from scratch
from_scratch <- dtw2vec(Q, c(C, y1, y2) , ws = WS)$normalized_distance
x$normalized_distance - from_scratch
plot(x)

# The plot shows alignments of high costs at the end 
# => attempt a decremtal step to find better partial matching
x <- decrement(x, direction = "C", refresh_dtw = TRUE)
x
plot(x)


#--- 2. example: First increment, then reverse increment
rw <- function(nn) cumsum(rnorm(nn))
Q <- rw(100)
C <- Q[11:90] + rnorm(80, 0, 0.1)
WS <- 40

# initial calculation
x <- initialize_plane(Q, C, ws = WS)
plot(x)

# incremental calculation for new observations that 
# are appened at the end of C
y1 <- Q[91:100] + rnorm(10, 0, 0.1)
x <- increment(x, newObs = y1)

# reverse the order of Q and C
x <- reverse(x)

# append new observations at the beginning: the new
# obervations must be in the same order as Q and C 
# => so newObs must be in reverse order, so y2 is 
# defined as Q from 10 to 6 (plus noise).
y2 <- Q[10:6] + rnorm(5, 0, 0.1)
x <- increment(x, newObs = y2)

# another incremental step in the reverse direction
y3 <- Q[5:1] + rnorm(5, 0, 0.1)
x <- increment(x, newObs = y3)

# compare with calculations from scratch, and plot x
from_scratch <- dtw2vec(rev(Q), rev(c(rev(y3), rev(y2), C, y1)),
                        ws = WS)$distance
x$distance - from_scratch
print(x)
plot(x)

}
}

\keyword{ ts }% __ONLY ONE__ keyword per line

