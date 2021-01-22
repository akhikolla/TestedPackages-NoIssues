\name{IntervalSurgeon-package}
\alias{IntervalSurgeon-package}
\alias{IntervalSurgeon}
\docType{package}
\title{\packageTitle{IntervalSurgeon}}
\description{\packageDescription{IntervalSurgeon}}
\details{
  The DESCRIPTION file: \packageDESCRIPTION{IntervalSurgeon}
  \packageIndices{IntervalSurgeon}

  IntervalSurgeon presents functions for manipulating integer-bounded sets of intervals.
  Sets of intervals are represented by two-column matrices, where inclusive start points are stored in the first column, and exclusive end points in the second.
  A central concept in the package is the `sections' of a set of intervals \code{x}: the non-overlapping, completely-covering set of intervals on the range of \code{x}, formed by making intervals between the consecutive sorted start/end points of the intevals in \code{x}.
  The function \code{sections} returns such a set of intervals given an input set.
  
}
\author{
  \packageAuthor{IntervalSurgeon}
  
  Maintainer: \packageMaintainer{IntervalSurgeon}
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{package}
\seealso{
  Optional links to other man pages
}
\examples{
  ## Optional simple examples of the most important functions
  ## Use \dontrun{} around code to be shown but not executed
}
