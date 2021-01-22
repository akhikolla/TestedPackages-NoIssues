\name{scaleMargins}
\alias{scaleMargins}
\title{
  Scale Rows and/or Columns of a Matrix (wordspace)
}
\description{

  This function provides a fast and memory-efficient way to scale the rows and/or columns
  of a dense or sparse matrix.  Each row is multiplied with the corresponding element of
  the vector \code{rows}, each column with the corresponding element of the vector \code{cols}.

}
\usage{

scaleMargins(M, rows=NULL, cols=NULL, duplicate=TRUE)

}
\arguments{

  \item{M}{a dense or sparse matrix in canonical format}

  \item{rows}{a numeric vector with length equal to the number of rows of M, or a single number.
    If missing or \code{NULL}, the rows of M are not rescaled.
  }

  \item{cols}{a numeric vector with length equal to the number of columns of M, or a single number.
    If missing or \code{NULL}, the columns of M are not rescaled.
  }

  \item{duplicate}{if \code{FALSE}, modify the matrix \code{M} in place.  Don't ever set this argument to \code{FALSE}.}
}
\details{
  
  If \code{M} is not in canonical format (dense numeric matrix or sparse matrix of class \code{dgCMatrix}),
  it will automatically be converted.  In this case, the precise behaviour of \code{duplicate=FALSE} is undefined.

  \code{duplicate=FALSE} is intended for internal use only.  
  
}
\value{

  The rescaled dense or sparse matrix.

}
\seealso{

  \code{\link{dsm.is.canonical}}

}
\examples{

M <- matrix(1, 5, 3)
scaleMargins(M, rows=1:5, cols=c(1, 10, 100))

}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}

