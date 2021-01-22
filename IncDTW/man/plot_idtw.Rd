\name{plot.idtw}
\alias{plot_idtw}
\alias{plot.idtw}
\alias{plot_planedtw}
\alias{plot.planedtw}
\alias{plotQC}
\alias{plotWarp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the results from Dynamic Time Warping
}
\description{
Plot function for objects of type \code{idtw}, the output of \code{dtw()} and \code{idtw()} respectively.}
\usage{
\method{plot}{idtw}(x, type = c("QC", "warp"), partial = NULL, selDim = 1, ...)
        
# an alias for plot_idtw
plot_idtw(x, type = c("QC", "warp"), partial = NULL, selDim = 1, ...)

\method{plot}{planedtw}(x, type = c("QC", "warp"), partial = NULL, selDim = 1, ...)

# an alias for plot_planedtw
plot_planedtw(x, type = c("QC", "warp"), partial = NULL, selDim = 1, ...)

plotQC(x, Q, C, partial = NULL, selDim = 1, ...)

plotWarp(x, Q, C, partial = NULL, selDim = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{output from \code{dtw(Q, C)}}
 \item{Q}{one dimensional numeric vector}
\item{C}{one dimensional numeric vector}
\item{type}{character, one of c('QC', 'warp')}
\item{partial}{list, the return value of \code{\link{dtw_partial}()}. Default = NULL, see \code{\link{dtw_partial}()} for details.}
\item{selDim}{integer, gives the column index of the multivariate time series (matrices) to be plotted. (default = 1) If Q and C are univariate time series (vectors) then selDim is neglected. }
\item{...}{Other arguments passed on to methods.}
}
\details{
The plot function visualizes the time warp and the alignment of the two time series. Also for partial alignments see \code{\link{dtw_partial}()}
}



\examples{

Q <- cumsum(rnorm(100))
C <- Q[11:100] + rnorm(90, 0, 0.5)
tmp <- dtw(Q = Q, C = C, ws = 15, return_wp = TRUE, return_QC = TRUE) 

plot(tmp, type = 'QC')
plotQC(tmp)
plot(tmp, type = 'warp')
plotWarp(tmp)


}
