\name{rsvd}
\alias{rsvd}
\title{
  Randomized Singular Value Decomposition (wordspace)
}
\description{

  An implementation of the randomized truncated SVD algorithm of Halko, Martinsson & Tropp (2009).

}
\usage{

rsvd(M, n, q = 2, oversampling = 2, transpose = FALSE, verbose = FALSE)

}
\arguments{

  \item{M}{
    a dense or sparse numeric matrix
  }

  \item{n}{
    an integer specifying the desired number of singular components. This argument must be specified and must satisfy \code{n <= min(nrow(M), ncol(M))}.
  }

  \item{q}{
    number of power iterations (Halko \emph{et al.} recommend \code{q=1} or \code{q=2})
  }

  \item{oversampling}{
    oversampling factor. The rSVD algorithm computes an approximate SVD factorization of rank \code{n * oversampling}, which is then truncated to the first \code{n} components.
  }

  \item{transpose}{
    if \code{TRUE}, apply the rSVD algorithm to the transpose \code{t(M)}, which may be more efficient depending on the dimensions of \code{M}
  }

  \item{verbose}{
    whether to display progress messages during execution
  }
}

\details{

  This implementation of randomized truncated SVD is based on the randomized PCA algorithm (Halko \emph{et al.} 2009, p. 9).  The discussion in Sec. 4 and 5 of the paper shows that the same algorithm applies to the case where the columns of A are not centered (Algorithm 4.3 + Algorithm 5.1).

}

\value{

  A list with components
  
  \item{u}{a matrix whose columns contain the first \code{n} left singular vectors of \code{M}}

  \item{v}{a matrix whose columns contain the first \code{n} right singular vectors of \code{M}}
  
  \item{d}{a vector containing the first \code{n} singular values of \code{M}}

}

\references{

  Halko, N., Martinsson, P. G., and Tropp, J. A. (2009).
  Finding structure with randomness: Stochastic algorithms for constructing approximate matrix decompositions. Technical Report 2009-05, ACM, California Institute of Technology.

}
\seealso{

  \code{\link{svd}}, \code{\link{dsm.projection}}, \code{\link{sparsesvd}}
  
}
%% \examples{
%% }
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
