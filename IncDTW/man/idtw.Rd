\name{idtw}
\alias{idtw}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Incremental DTW
}
\description{
Update the DTW distance, cost matrices and direction matrices including the warping path for new observations of two time series.
}
\usage{
idtw(Q, C, newObs, gcm, dm, 
      dist_method = c("norm1", "norm2", "norm2_square"),
      step_pattern = c("symmetric2", "symmetric1"),
      diffM = NULL, ws = NULL, 
      return_cm = FALSE,
      return_diffM = FALSE,
      return_wp = FALSE,
      return_diffp = FALSE,
      return_QC = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Q}{numeric vector, or matrix (see also \code{\link{dtw}})}
\item{C}{numeric vector, or matrix}
\item{newObs}{vector or matrix of new observations to be appended to C}
\item{gcm}{global cost matrix, output from \code{dtw(Q, C, ...)}}
\item{dm}{direction matrix, output from \code{dtw(Q, C, ...)}}
\item{dist_method}{character, describes the method of distance measure. See also \code{\link{dtw}}.}
\item{step_pattern}{character, describes the step pattern. See also \code{\link{dtw}}.}
\item{diffM}{differences matrix, output from \code{dtw(Q, C, ...)}. This matrix is an optional input parameter (default = NULL) that is necessary to return the path of differences. Only for univariate time series Q and C.}
\item{ws}{integer, describes the window size for the sakoe chiba window. If NULL, then no window is applied. (default = NULL)
}
\item{return_cm}{logical, if TRUE then the Matrix of costs (the absolute value) is returned. (default = FALSE)
}
\item{return_diffM}{logical, if TRUE then the Matrix of differences (not the absolute value) is returned. (default = FALSE)
}
\item{return_wp}{
logical, if TRUE then the warping path is returned. (default = FALSE) If return_diffp == TRUE, then return_wp is set to TRUE as well.
}
\item{return_diffp}{
logical, if TRUE then the path of differences (not the absolute value) is returned. (default = FALSE)
}
\item{return_QC}{
logical, if TRUE then the input vectors Q and C are appended to the returned list. This is useful for the \code{\link{plot.idtw}} function.  (default = FALSE)}
}
\details{
The dynamic time warping distance is the element in the last row and last column of the global cost matrix.
}
\value{
\item{distance}{the DTW distance, that is the element of the last row and last column of gcm}
\item{gcm}{global cost matrix}
\item{dm}{direction matrix (3=up, 1=diagonal, 2=left)}
\item{wp}{warping path}
\item{ii}{indices of Q of the optimal path}
\item{jj}{indices of C of the optimal path}
\item{cm}{Matrix of costs}
\item{diffM}{Matrix of differences}
\item{diffp}{path of differences}
\item{Q}{input Q}
\item{C}{input C}
\item{normalized_distance}{the normalized DTW distance, see also \code{link{dtw}}}
}
\references{
Sakoe, H.; Chiba, S., Dynamic programming algorithm optimization for spoken word recognition, Acoustics, Speech, and Signal Processing [see also IEEE Transactions on Signal Processing], IEEE Transactions on , vol.26, no.1, pp. 43-49, Feb 1978.
http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055
}


\examples{
#--- Compare the incremental calculation with the basic
# calculation from scratch.
Q <- cumsum(rnorm(100))
C <- Q[11:100] + rnorm(90, 0, 0.5)
newObs <-  c(2, 3)# new observation
base <- dtw(Q = Q, C = C, ws = 15, return_diffM = TRUE) 
base

# recalculation from scratch with new observations
result0 <- dtw(Q = Q, C = c(C, newObs), ws = 15,  return_diffM = TRUE)

# the incremental step with new observations
result1 <- idtw(Q, C, ws = 15, newO = newObs, gcm = base$gcm, 
                dm = base$dm, diffM = base$diffM, return_diffp = TRUE,  
                return_diffM = TRUE, return_QC = TRUE) 
print(result1, digits = 2)
plot(result1)


#--- Compare the incremental calculation with external calculated 
# costMatrix cm_add with the basic calculation from scratch.
cm_add <- cm(Q, newObs)
result2 <- idtw(Q = cm_add, C = "cm_add", ws = 15, newO = newObs, 
                gcm = base$gcm, dm = base$dm) 

c(result0$distance, result1$distance, result2$distance)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }% use one of  RShowDoc("KEYWORDS")
\keyword{ ts }% __ONLY ONE__ keyword per line
