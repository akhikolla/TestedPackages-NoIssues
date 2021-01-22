\name{dsm}
\alias{dsm}
\alias{dsm-object}
\title{
  Create DSM Object Representing a Distributional Semantic Model (wordspace)
}
\description{

  This is the constructor function for \code{dsm} objects representing distributional semantic models,
  i.e. a co-occurrence matrix together with additional information on target terms (rows) and features (columns).
  A new DSM can be initialised with a dense or sparse co-occurrence matrix, or with a triplet representation of a sparse matrix.

}
\usage{

dsm(M = NULL, target = NULL, feature = NULL, score = NULL,
    rowinfo = NULL, colinfo = NULL, N = NA,
    globals = list(), raw.freq = FALSE, sort = FALSE, verbose = FALSE)

}
\arguments{

  \item{M}{
    a dense or sparse co-occurrence matrix. A sparse matrix must be a subclass of \code{\link[=sparseMatrix-class]{sparseMatrix}} from the \code{Matrix} package. See "Details" below.
  }

  \item{target}{
    a character vector of target terms (see "Details" below)
  }

  \item{feature}{
    a character vector of feature terms (see "Details" below)
  }

  \item{score}{
    a numeric vector of co-occurrence frequencies or weighted/transformed scores (see "Details" below)
  }

  \item{rowinfo}{
    a data frame containing information about the rows of the co-occurrence matrix, corresponding to target terms.  The data frame must include a column \code{term} with the target term labels.  If unspecified, a minimal \code{rowinfo} table is compiled automatically (see "Details" below).
  }

  \item{colinfo}{
    a data frame containing information about the columns of the co-occurrence matrix, corresponding to feature terms.  The data frame must include a column \code{term} with the feature term labels.  If unspecified, a minimal \code{colinfo} table is compiled automatically (see "Details" below).
  }

  \item{N}{
    a single numeric value specifying the effective sample size of the co-occurrence matrix.  This value may be determined automatically if \code{raw.freq=TRUE}.
  }

  \item{globals}{
    a list of global variables, which are included in the \code{globals} field of the DSM object.  May contain an entry for the sample size \eqn{N}, which can be overridden by an explicitly specified value in the argument \code{N}.
  }
  
  \item{raw.freq}{
    if \code{TRUE}, entries of the co-occurrence matrix are interpreted as raw frequency counts. By default, it is assumed that some weighting/transformation has already been applied.
  }

  \item{sort}{
    if \code{TRUE}, sort rows and columns of a co-occurrence matrix specified in triplet form alphabetically. If the matrix is given directly (in argument \code{M}), rows and columns are never reordered.
  }
  
  \item{verbose}{
    if \code{TRUE}, a few progress and information messages are shown
  }

}
\details{
  
  The co-occurrence matrix forming the core of the distributional semantic model (DSM) can be specified in two different ways:

  \enumerate{
    \item{
      As a dense or sparse matrix in argument \code{M}.  A sparse matrix must be a subclass of \code{\link[=dMatrix-class]{dMatrix}} (from the \code{Matrix} package) and is automatically converted to the canonical storage mode used by the \code{wordspace} package.  Row and column labels may be specified with arguments \code{target} and \code{feature}, which must be character vectors of suitable length; otherwise \code{dimnames(M)} are used.
    }
    \item{
      As a triplet representation in arguments \code{target} (row label), \code{feature} (column label) and \code{score} (co-occurrence frequency or pre-computed score).  The three arguments must be vectors of the same length; each set of corresponding elements specifies a non-zero cell of the co-occurrence matrix.  If multiple entries for the same cell are given, their frequency or score values are added up.
    }
  }

  The optional arguments \code{rowinfo} and \code{colinfo} are data frames with additional information about target and feature terms.  If they are specified, they must contain a column \code{$term} matching the row or column labels of the co-occurrence matrix.  Marginal frequencies and nonzero or document counts can be given in columns \code{$f} and \code{$nnzero}; any further columns are interpreted as meta-information on the target or feature terms.  The rows of each data frame are automatically reordered to match the rows or columns of the co-occurrence matrix.  Target or feature terms that do not appear in the co-occurrence matrix are silently discarded.

  Counts of nonzero cells for each row and column are computed automatically, unless they are already present in the \code{rowinfo} and \code{colinfo} data frames.  If the co-occurrence matrix contains raw frequency values, marginal frequencies for the target and feature terms are also computed automatically unless given in \code{rowinfo} and \code{colinfo}; the same holds for the effective sample size \code{N}.

  If \code{raw.freq=TRUE}, all matrix entries must be non-negative; fractional frequency counts are allowed, however.
  
}
\value{

  An object of class \code{dsm}, a list with the following components:
  
  \item{M}{
    A co-occurrence matrix of raw frequency counts in canonical format (see \code{\link{dsm.canonical.matrix}}).
  }
  
  \item{S}{
    A weighted and transformed co-occurrence matrix ("score" matrix) in canonical format (see \code{\link{dsm.canonical.matrix}}).
    Either \code{M} or \code{S} or both may be present.  The object returned by \code{dsm()} will include \code{M} if \code{raw.freq=TRUE} and \code{S} otherwise.
  }

  \item{rows}{
    A data frame with information about the target terms, corresponding to the rows of the co-occurrence matrix.  The data frame usually has at least three columns:
    \describe{
      \item{\code{rows$term}}{the target term = row label}
      \item{\code{rows$f}}{marginal frequency of the target term; must be present if the DSM object contains a raw co-occurrence matrix \code{M}}
      \item{\code{rows$nnzero}}{number of nonzero entries in the corresponding row of the co-occurrence matrix}
    }
    Further columns may provide additional information.
  }
  
  \item{cols}{
    A data frame with information about the feature terms, corresponding to the columns of the co-occurrence matrix, in the same format as \code{rows}.
  }

  \item{globals}{
    A list of global variables.  The following variables have a special meaning:
    \describe{
      \item{\code{globals$N}}{effective sample size of the underlying corpus; may be \code{NA} if raw co-occurrence counts are not available}
      \item{\code{globals$locked}}{if \code{TRUE}, the marginal frequencies are no longer valid due to a \code{merge}, \code{rbind} or \code{cbind} operation; in this case, association scores cannot be computed from the co-occurrence frequencies \code{M}}
    }
  }

}
\seealso{

  See \code{\link{dsm.canonical.matrix}} for a description of the canonical matrix formats.  DSM objects are usually loaded directly from a disk file in UCS (\code{\link{read.dsm.ucs}}) or triplet (\code{\link{read.dsm.triplet}}) format.

}
\examples{

MyDSM <- dsm(
  target =  c("boat", "boat", "cat",  "dog", "dog"),
  feature = c("buy",  "use",  "feed", "buy", "feed"),
  score =   c(1,      3,      2,      1,     1),
  raw.freq = TRUE
)

print(MyDSM)   # 3 x 3 matrix with 5 out of 9 nonzero cells
print(MyDSM$M) # the actual co-occurrence matrix

print(MyDSM$rows) # row information
print(MyDSM$cols) # column information

}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{classes}
