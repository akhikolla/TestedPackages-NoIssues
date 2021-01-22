\name{dtw2vec}
\alias{dtw2vec}
\alias{dtw2vec_cm}
\alias{dtw2vec_univ}
\alias{dtw2vec_multiv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast vector-based Dynamic Time Warping
}
\description{
Calculates the Dynamic Time Warping distance by hand of a vector-based implementation and is much faster than the traditional method \code{dtw()}. Also allows early abandoning and sakoe chiba warping window, both for univariate and multivariate time series.
}
\usage{
dtw2vec(Q, C, dist_method = c("norm1", "norm2", "norm2_square"),
         step_pattern = c("symmetric2", "symmetric1"),
         ws = NULL, threshold = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Q}{Either Q is (a) a time series (vector or matrix for multivariate time series) or (b) Q is a cost matrix, so a matrix storing the local distances of the time series Q and C. If Q and C are matrices, they need to have the same number of columns. If Q is a cost matrix, C needs to be equal the character string "cm".}
\item{C}{time series as vector or matrix, or for case (b) C equals "cm"}
\item{dist_method}{character, describes the method of distance measure. See also \code{\link{dtw}}.  If Q is a cost matrix, the dist_method parameter is not necessary.}
\item{step_pattern}{character, describes the step pattern. See also \code{\link{dtw}}.}
\item{ws}{
integer, describes the window size for the sakoe chiba window. If NULL, then no window is applied. (default = NULL)
}
\item{threshold}{numeric, the threshold for early abandoning. In the calculation of the global cost matrix a possible path stops as soon as the threshold is reached. Facilitates faster calculations in case of low threshold. The threshold relates to the non-normalized distance measure. (default = NULL, no early abandoning)
}
}

\details{
no matrices are allocated, no matrices are returned
}
\value{
\item{distance}{the DTW distance}
\item{normalized_distance}{the normalized DTW distance, see also \code{\link{dtw}}}
}
\references{
Sakoe, H.; Chiba, S., Dynamic programming algorithm optimization for spoken word recognition, Acoustics, Speech, and Signal Processing [see also IEEE Transactions on Signal Processing], IEEE Transactions on , vol.26, no.1, pp. 43-49, Feb 1978. http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163055
}


\examples{
Q <- cumsum(rnorm(100))
C <- Q[11:100] + rnorm(90, 0, 0.5)
dtw2vec(Q = Q, C = C)
dtw2vec(Q = Q, C = C, ws = 30)
dtw2vec(Q = Q, C = C, threshold = 100)
dtw2vec(Q = Q, C = C, ws = 30, threshold = 100)

cm0 <- cm(Q, C)
dtw2vec(Q = cm0, C = "cm", ws = 30, threshold = 100)
}
\keyword{ classif }% use one of  RShowDoc("KEYWORDS")
\keyword{ cluster }% use one of  RShowDoc("KEYWORDS")
\keyword{ ts }% __ONLY ONE__ keyword per line

