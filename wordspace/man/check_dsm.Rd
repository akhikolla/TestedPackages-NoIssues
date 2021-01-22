\name{check.dsm}
\alias{check.dsm}
\title{
  Validate Internal Structure of DSM Object (wordspace)
}
\description{

  Validate the internal structure of a DSM object and return a list with information about the object.

}
\usage{

check.dsm(model, validate = FALSE, nonneg.check = FALSE)

}
\arguments{

  \item{model}{
    an object of class \code{dsm}
  }
  
  \item{validate}{
    carry out extended validation of internal consistency? (may be expensive)
  }

  \item{nonneg.check}{
    if \code{TRUE}, check the co-occurrence (\eqn{M}) and/or score (\eqn{S}) matrix for non-negativity (may be expensive)
  }
  
}
%% \details{
%% }
\value{

  Aborts with error message if any inconsistency is detected.
  Otherwise a list with the following items is returned:

  \item{nrow}{
    number of rows (target terms) of the DSM
  }
  \item{ncol}{
    number of columns (features) of the DSM
  }
  \item{N}{
    sample size of the underlying data set (may be \code{NA})
  }
  \item{M$ok}{
    whether co-occurrence frequency matrix \eqn{M} is available
  }
  \item{M$sparse}{
    whether \eqn{M} is sparse or dense (only present if \code{M$ok})
  } 
  \item{M$canonical}{
    whether \eqn{M} is in canonical DSM format (only present if \code{M$ok})
  }
  \item{M$nonneg}{
    whether \eqn{M} is non-negative (only present if \code{M$ok}, and may be \code{NA} unless \code{nonneg.check=TRUE was specified})
  }
  \item{S$ok}{
    whether score matrix \eqn{S} is available
  }
  \item{S$sparse}{
    whether \eqn{S} is sparse or dense (only present if \code{S$ok})
  } 
  \item{S$canonical}{
    whether \eqn{S} is in canonical DSM format (only present if \code{S$ok})
  }
  \item{S$nonneg}{
    whether \eqn{S} is non-negative (only present if \code{S$ok}, and may be \code{NA} unless \code{nonneg.check=TRUE was specified})
  }
  \item{locked}{
    \code{TRUE} if matrix combines data with inconsistent row or column marginals (in this case, association scores cannot be computed any more)
  }
}
\seealso{
  \code{\link{dsm}}, \code{\link{print.dsm}}
}
\examples{

check.dsm(DSM_TermTerm)

}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
