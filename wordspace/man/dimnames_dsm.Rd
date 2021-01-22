\name{dimnames.dsm}
\alias{dimnames.dsm}
\alias{dimnames<-.dsm}
\title{
  Dimnames of a DSM Object (wordspace)
}
\description{

  Retrieve or set the dimnames of the co-occurrence and/or score matrix
  represented by a DSM object.

}
\usage{

\method{dimnames}{dsm}(x)
\method{dimnames}{dsm}(x) <- value

}
\arguments{

  \item{x}{
    an object of class \code{dsm}
  }
  
  \item{value}{
    a list of two character vectors with new row and column names for \code{x}. Both vectors must have appropriate length and may not be \code{NULL}.
  }

}
\details{

  This method automatically checks that the row and column names of the co-occurrence and/or score matrix are consistent with the target terms and features listed in the row/column information tables.
  

}
\value{

  The \code{dimnames()} of a DSM object are always a list of length 2, consisting of two character vectors with row and column labels, respectively.

}
\seealso{

  \code{\link{dsm}}, \code{\link{check.dsm}}, \code{\link{dim.dsm}}
  
}
\examples{
rownames(DSM_TermContext)
colnames(DSM_TermContext)

tmp <- DSM_TermContext
rownames(tmp)[3] <- "pet"
head(tmp, 4, 6)
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
