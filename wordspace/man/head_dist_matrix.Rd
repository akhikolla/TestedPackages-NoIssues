\name{head.dist.matrix}
\alias{head.dist.matrix}
\title{
  Return the Top Left Corner of a Distance Matrix (wordspace)
}
\description{

  Returns the first \code{n} rows and first \code{k} columns of a distance matrix returned by the \code{\link{dist.matrix}} function.

}
\usage{

\method{head}{dist.matrix}(x, n = 6L, k = n, \dots)

}
\arguments{

  \item{x}{
    an distance matrix of class \code{dist.matrix}
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

  A numeric matrix with \code{n} rows and \code{k} columns.

  \bold{Note:} this matrix is no longer marked as an object of class \code{dist.matrix} and thus prints nicely without attributes.

}
\seealso{

  \code{\link{head}} for the generic method.
  
}
\examples{

dm <- dist.matrix(DSM_Vectors[1:100, ])
print(head(dm, 8, 5), digits=3)

}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
