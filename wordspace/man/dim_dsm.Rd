\name{dim.dsm}
\alias{dim.dsm}
\title{
  Dimensions of a DSM Object (wordspace)
}
\description{

  Retrieve the dimensions of the co-occurrence and/or score matrix 
  represented by a DSM object.

}
\usage{

\method{dim}{dsm}(x)

}
\arguments{

  \item{x}{
    an object of class \code{dsm}
  }

}
\details{

  Note that an assignment form (\code{dim<-}) for modifying dimensions is not provided.

}
\value{

  An integer vector of length 2, specifying the number of rows and the number of columns of the DSM matrix.
  
}
\seealso{

  \code{\link{dsm}}, \code{\link{check.dsm}}, \code{\link{print.dsm}}, \code{\link{dimnames.dsm}}
  
}
\examples{
dim(DSM_TermTerm)
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
