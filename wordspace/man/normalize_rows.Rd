\name{normalize.rows}
\alias{normalize.rows}
\alias{normalize.cols}
\title{
  Normalize Rows or Columns of Matrix to Unit Length (wordspace)
}
\description{

  Efficiently normalize the row or column vectors of a dense or sparse matrix to unit length.

}
\usage{

normalize.rows(M, method = "euclidean", p = 2, \dots, 
               tol = 1e-6, inplace = FALSE)

normalize.cols(M, method = "euclidean", p = 2, \dots,
               tol = 1e-6, inplace = FALSE)

}
\arguments{

  \item{M}{a dense or sparse numeric matrix}

  \item{method}{norm to be computed, see \code{\link{rowNorms}}}

  \item{p}{exponent of Minkowski p-norm in the range \eqn{0 < p \le \infty}. Note that normalization is not possible for very small values of \eqn{p}.}

  \item{\dots}{any further arguments are passed to \code{\link{rowNorms}} (or \code{\link{colNorms}})}

  \item{tol}{row/column vectors with norm below \code{tol} are assumed to be all zeroes and cannot be normalized (see \dQuote{Details} below)}
  
  \item{inplace}{if \code{TRUE}, modify the matrix \code{M} in place.  Don't ever set this argument to \code{TRUE}.}
}
\details{
  These functions return a matrix with row (or column) vectors rescaled to a length of 1 according to the selected norm.
  
  All-zero vectors (with \eqn{\|0\| = 0}{|0| = 0}) cannot be normalized.  In order to avoid scaling up rounding errors, rows (or columns) with \eqn{\|x\| < tol}{|x| < tol} are explicitly set to 0 (and thus not normalized). Since a suitable threshold for rounding errors depends on the scaling behaviour of the selected norm and the provenance of \eqn{M}, it is advisable to set \code{tol} explicitly to an appropriate value. Pass \code{tol = 0} to normalize all nonzero vectors.

  The generalized Minkowski norm with \eqn{p < 1} is not homogeneous but can still be normalized.  This is numerically unstable for very small values of \eqn{p}, which will be rejected with an error message.  The Hamming length (\eqn{p = 0}) cannot be normalized at all.  See \code{\link{rowNorms}} for more information.
}
\value{

  A row-normalized (or column-normalized) matrix with the same dimensions as \eqn{M}.

}
\seealso{
  See \code{\link{rowNorms}} for details on available norms and their parameters.
}
%% \examples{
%% }
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
