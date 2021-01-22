\name{head.dsm}
\alias{head.dsm}
\title{
  Return the Top Left Corner of a DSM Matrix (wordspace)
}
\description{

  Returns the first \code{n} rows and first \code{k} columns of the co-occurrence matrix stored in a \code{dsm} object.
  If a scored matrix is available, it is automatically used; otherwise the raw frequencies are shown.

}
\usage{

\method{head}{dsm}(x, n = 6L, k = n, \dots)

}
\arguments{

  \item{x}{
    an object of class \code{dsm}
  }

  \item{n}{
    a single integer specifying the number of rows to extract
  }

  \item{k}{
    a single integer specifying the number of columns to extract (default: same as number of rows)
  }

  \item{\dots}{
    all other arguments are silently ignored
  }
  
}
\details{

  Note that in contrast to other \code{\link{head}} methods, negative values of \code{n} (and \code{k}) are not supported.  There is also currently no corresponding \code{\link{tail}} method.

}
\value{

  A dense or sparse co-occurrence matrix with \code{n} rows and \code{k} columns.

}
\seealso{

  \code{\link{head}} for the generic method.
  
}
\examples{

head(DSM_TermTerm, Inf, Inf) # show full co-occurrence matrix

head(DSM_TermTerm, 3, 4)

}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
