\name{IncDTW-package}
\alias{IncDTW-package}
%\alias{idtw}
% https://colinfay.me/writing-r-extensions/writing-r-documentation-files.html
\docType{package}
\title{
  Incremental Dynamic Time Warping
}
\description{
The Dynamic Time Warping (DTW) distance for time series allows non-linear alignments of time series to match  similar patterns in time series of different lengths and or different speeds. Beside the traditional implementation of the DTW algorithm, the specialties of this package are, (1) the incremental calculation, which is specifically useful for life data streams due to computationally efficiency, (2) the vector based implementation of the traditional DTW algorithm which is faster because no matrices are allocated and is especially useful for computing distance matrices of pairwise DTW distances for many time series  and (3) the combination of incremental and vector-based calculation.
}

\details{
Main features:
\itemize{
\item{
Incremental Calculation, \code{\link{idtw}}, \code{\link{idtw2vec}} and \code{\link{increment}}
}
\item{
Detect k-nearest subsequences in longer time series, \code{\link{rundtw}}
}
\item{
Matrix-based \code{\link{dtw}} and Vector-based \code{\link{dtw2vec}} implementation of the DTW algorithm
}
\item{
Sakoe Chiba warping window
}
\item{
Early abandoning and lower bounding
}
\item{
Support for multivariate time series
}
\item{
Fast calculation of a distance matrix of pairwise DTW distances for clustering or classification of many multivariate time series, \code{\link{dtw_dismat}}
}
\item{
Aggregate cluster members with \code{\link{dba}} or get the centroid with \code{\link{centroid}}
}
\item{
C++ in the heart thanks to Rcpp
}
}

}

\author{
Maximilian Leodolter

Maintainer: Maximilian Leodolter <maximilian.leodolter@gmail.com>
}
\references{
  Dynamic programming algorithm optimization for spoken word recognition by Sakoe and Chiba published in 1978 (DOI:10.1109/TASSP.1978.1163055)
}
\keyword{ DTW }
\seealso{
\url{https://ieeexplore.ieee.org/document/1163055/}

\url{https://en.wikipedia.org/wiki/Dynamic_time_warping}

\url{https://github.com/maxar/IncDTW}
}
