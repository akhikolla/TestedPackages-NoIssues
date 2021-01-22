\name{signcount}
\alias{signcount}
\title{
  Efficiently Count Positive, Negative and Zero Values (wordspace)
}
\description{

  This function counts the number of positive, negative and zero elements in a numeric vector or matrix,
  including some types of dense and sparse representations in the \code{\link{Matrix}} package.
  
  It can be used to test for non-negativity and compute nonzero counts efficiently, without any memory overhead.
  
}
\usage{

signcount(x, what = c("counts", "nonneg", "nnzero"))

}
\arguments{

  \item{x}{
    a numeric vector or array, or a numeric \code{Matrix} (supported formats are \code{\link[=dgeMatrix-class]{dgeMatrix}},  \code{\link[=dgCMatrix-class]{dgCMatrix}} and  \code{\link[=dgRMatrix-class]{dgRMatrix}})
  }

  \item{what}{
    whether to return the counts of positive, negative and zero elements (\code{counts}), the number of nonzero elements (\code{nnzero}), or to test for non-negativity (\code{nonneg}); see \sQuote{Value} below
  }

}
\details{
  
  \code{x} must not contain any undefined values; \code{signcount} does not check whether this is actually the case.
  
}
\value{
  
  \describe{
    \item{\code{what="counts"}}{
      A labelled numeric vector of length 3 with the counts of positive (\code{pos}), zero (\code{zero}) and negative (\code{neg}) values.
    } 
    \item{\code{what="nonneg"}}{
      A single logical value: \code{TRUE} if \code{x} is non-negative, \code{FALSE} otherwise.
    }
    \item{\code{what="nonzero"}}{
      A single numeric value, the total number of nonzero elements in \code{x}.
    }
  }
}
\examples{
signcount(DSM_TermTermMatrix)     # dense matrix
signcount(DSM_TermContextMatrix)  # sparse dgCMatrix
signcount(DSM_TermContextMatrix, "nonneg") # confirm that it is non-negative
}
\keyword{utilities}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
