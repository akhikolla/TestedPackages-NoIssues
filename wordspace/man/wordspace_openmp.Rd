\name{wordspace.openmp}
\alias{wordspace.openmp}
\title{
  Control multi-core processing in wordspace functions (wordspace)
}
\description{

  Control whether multi-core processing is used by wordspace functions (if available)
  and how many threads are run in parallel.  See "Details" below for more information.
  
}
\usage{

wordspace.openmp(threads = NULL)

}
\arguments{

  \item{threads}{
    if specified, number of parallel threads to be used for multi-core processing
  }

}
\details{

  The \code{wordspace} package has experimental support for multi-core processing using OpenMP on some platforms.
  So far, only the \code{dist.matrix} function uses multi-core processing (for all distance measures except \code{cosine}).
  
  Even where supported, OpenMP is not enabled by default and has to be activated explicitly with \code{wordspace.openmp(threads=N)}, where \code{N} is
  the number of parallel threads to be used.
  
  Call \code{wordspace.openmp()} without arguments to find out whether OpenMP is supported on your platform
  and obtain information about the maximum number of threads available as well as the current setting.

  Note that multi-threading of other R packages and functions (such as optimised matrix algebra in the BLAS library) is never affected by this function.
  
}
\value{

  If \code{threads} is unspecified or \code{NULL}, a data frame with a single row and the following information is returned:

  \item{available}{\code{TRUE} if OpenMP multi-core support is available}
  \item{max}{maximum number of threads that can be used (0 if not available)}
  \item{threads}{currently selected number of threads (defaults to 1 if not available)}
}
%% \references{
%% }
\seealso{
  \code{\link{dist.matrix}}
}
\examples{
wordspace.openmp()
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
