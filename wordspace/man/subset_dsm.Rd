\name{subset.dsm}
\alias{subset.dsm}
\title{
  Subsetting Distributional Semantic Models (wordspace)
}
\description{

  Filter the rows and/or columns of a DSM object according to user-specified conditions.

}
\usage{

\method{subset}{dsm}(x, subset, select, recursive = FALSE, drop.zeroes = FALSE,
       matrix.only = FALSE, envir = parent.frame(), run.gc = FALSE, \dots)

}
\arguments{

  \item{x}{
    an object of class \code{dsm}
  }

  \item{subset}{
    Boolean expression or index vector selecting a subset of the rows; the expression can use variables \code{term} and \code{f} to access target terms and their marginal frequencies, \code{nnzero} for the number of nonzero elements in each row, further optional variables from the row information table, as well as global variables such as the sample size \code{N}
  }

  \item{select}{
    Boolean expression or index vector selecting a subset of the columns; the expression can use variables \code{term} and \code{f} to access feature terms and their marginal frequencies, \code{nnzero} for the number of nonzero elements in each column, further optional variables from the column information table, as well as global variables such as the sample size \code{N}
  }

  \item{recursive}{if \code{TRUE} and both \code{subset} and \code{select} conditions are specified, the \code{subset} is applied repeatedly until the DSM no longer changes.
    This is typically needed if conditions on nonzero counts or row/column norms are specified, which may be affected by the subsetting procedure.
  }
  
  \item{drop.zeroes}{
    if \code{TRUE}, all rows and columns without any nonzero entries after subsetting are removed from the model
    (nonzero counts are based on the score matrix \eqn{S} if available, raw cooccurrence frequencies \eqn{M} otherwise)
  }

  \item{matrix.only}{
    if \code{TRUE}, return only the selected subset of the score matrix \eqn{S} (if available) or frequency matrix \eqn{M}, not a full DSM object.  This may conserve a substantial amount of memory when processing very large DSMs.
  }

  \item{envir}{
    environment in which the \code{subset} and \code{select} conditions are evaluated.  Defaults to the context of the function call, so all variables visible there can be used in the expressions.
  }

  \item{run.gc}{whether to run the garbage collector after each iteration of a recursive subset (\code{recursive=TRUE}) in order to keep memory overhead as low as possible. This option should only be specified if memory is very tight, since garbage collector runs can be expensive (e.g. when there are many distinct strings in the workspace).}

  \item{\dots}{
    any further arguments are silently ignored
  }
}
%% \details{
%% }
\value{

  An object of class \code{dsm} containing the specified subset of the model \code{x}.

  If necessary, counts of nonzero elements for each row and/or column are updated automatically.

}
\seealso{

  \code{\link{dsm}}

}
\examples{

print(DSM_TermContext$M)
model <- DSM_TermContext

subset(model, nchar(term) <= 4)$M     # short target terms
subset(model, select=(nnzero <= 3))$M # columns with <= 3 nonzero cells

subset(model, nchar(term) <= 4, nnzero <= 3)$M # combine both conditions

subset(model, nchar(term) <= 4, nnzero >= 2)$M # still three columns with nnzero < 2
subset(model, nchar(term) <= 4, nnzero >= 2, recursive=TRUE)$M

}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
