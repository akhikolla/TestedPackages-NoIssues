\name{dsm.canonical.matrix}
\alias{dsm.is.canonical}
\alias{dsm.canonical.matrix}
\title{
  Canonical Formats for a DSM Co-occurrence Matrix (wordspace)
}
\description{

  Test whether a co-occurrence matrix is represented in a DSM canonical
  format, or convert matrix to canonical format.
  
}
\usage{

dsm.is.canonical(x, nonneg.check = FALSE)

dsm.canonical.matrix(x, triplet = FALSE, annotate = FALSE, nonneg.check = FALSE)

}
\arguments{

  \item{x}{
    a dense or sparse DSM co-occurrence matrix
  }

  \item{nonneg.check}{
    if \code{TRUE}, check whether all elements of the matrix are non-negative
  }

  \item{triplet}{
    if \code{TRUE} and if \code{x} is sparse, return a matrix in triplet format (class \code{\link[=dgTMatrix-class]{dgTMatrix}}) rather than in column-compressed format (class \code{\link[=dgCMatrix-class]{dgCMatrix}}).  Note that this is \emph{not} a canonical DSM format.
  }

  \item{annotate}{
    if \code{TRUE}, annotate \code{x} with attributes \code{sparse} and \code{nonneg}, indicating whether the matrix is in sparse representation and non-negative, respectively.  Non-negativity is only checked if \code{nonneg.check=TRUE}; otherwise an existing attribute will be passed through without validation.
  }

}
\details{
  
  Note that conversion into canonical format may result in unnecessary copying of \code{x}, especially if \code{annotate=TRUE}.
  For optimal performance, set \code{annotate=FALSE} whenever possible and do not call \code{dsm.canonical.matrix()} as a no-op.

  Instead of \preformatted{    M <- dsm.canonical.matrix(M, annotate=TRUE, nonneg=TRUE)}
  use \preformatted{    M.flags <- dsm.is.canonical(M, nonneg=FALSE)
    if (!M.flags$canonical) M <- dsm.canonical.matrix(M)
    M.flags <- dsm.is.canonical(M, nonneg=TRUE)}

  If \code{nonneg.check=FALSE} and \code{x} has an attribute \code{nonneg}, its value is accepted without validation.
  
  Checking non-negativity can be expensive and create substantial memory overhead.  It is guaranteed to be efficient for a matrix in canonical format.
}
\value{

  \code{dsm.is.canonical()} returns a data frame containing a single row with the following items:

  \item{sparse}{whether \code{x} is a sparse (\code{TRUE}) or dense (\code{TRUE}) matrix}
  \item{canonical}{whether \code{x} is in canonical format}
  \item{nonneg}{whether all cells of \code{x} are non-negative; may be \code{NA} if \code{nonneg.check=FALSE}}
  
  \code{dsm.canonical.matrix()} returns a matrix in canonical DSM format, i.e.
  \itemize{
    \item{of class \code{\link{matrix}} for a dense matrix (even if \code{x} is a \code{\link[=denseMatrix-class]{denseMatrix}} object);}
    \item{of class \code{\link[=dgCMatrix-class]{dgCMatrix}} for a sparse matrix.}
  }
  If \code{triplet=TRUE} and \code{x} is sparse, it returns a matrix of class \code{\link[=dgTMatrix-class]{dgTMatrix}}, which is \emph{not} a canonical format.

  If \code{annotate=TRUE}, the returned matrix has attributes \code{sparse} and \code{nonneg} (possibly \code{NA}).
}
%% \seealso{
%% }
%% \examples{
%% }
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
