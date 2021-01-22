\name{as.matrix.dsm}
\alias{as.matrix.dsm}
\title{
  Extract Matrix from DSM Object (wordspace)
}
\description{

  Extract the co-occurrence or score matrix from a DSM object.

}
\usage{

\method{as.matrix}{dsm}(x, what = c("auto", "M", "S"), ...)

}
\arguments{

  \item{x}{
    an object of class \code{dsm}
  }
  
  \item{what}{
    whether to extract the raw co-occurrence matrix (\code{M}) or the score matrix (\code{S}). The default option \code{auto} prefers the score matrix if both are available.
  }

  \item{\dots}{
    any additional arguments are ignored
  }

}
\details{

  This function ensures that the row and column names of the matrix are consistent with the row/column information tables of the DSM.
  For faster access to the matrix, simply use \code{x$M} or \code{x$S} directly.

}
\value{

  Either the raw co-occurrence matrix or the score matrix of the DSM \code{x}.
  
  Note that unlike other \code{as.matrix} methods, a sparse matrix in canonical DSM format may be returned.

}
\seealso{

  \code{\link{dsm}}, \code{\link{dim.dsm}}, \code{\link{dimnames.dsm}}, \code{\link{dsm.is.canonical}}
  
}
\examples{
as.matrix(DSM_TermTerm)
as.matrix(DSM_TermContext)
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
