\name{nearest.neighbours}
\alias{nearest.neighbours}
\alias{nearest.neighbors}
\title{
  Find Nearest Neighbours in DSM Space (wordspace)
}
\description{

  Find the nearest neighbours of a term vector in a DSM, given either as a scored cooccurrence matrix or a pre-computed distance matrix.  The target term can be selected by name (in which case the cooccurrence or distance matrix must be labelled appropriately) or specified as a vector (if the DSM is given as a matrix).

}
\usage{

nearest.neighbours(M, term, n = 10, M2 = NULL, byrow = TRUE,
                   drop = TRUE, skip.missing = FALSE, dist.matrix = FALSE,
                   \dots, batchsize=50e6, verbose=FALSE)

}
\arguments{

  \item{M}{
    either a dense or sparse matrix representing a scored DSM (or an object of class \code{dsm}), or a pre-computed distance matrix returned by \code{dist.matrix} (as an object of class \code{dist.matrix}).  Note that the compact representation produced by the \code{dist} function (class \code{dist}) is not accepted.
  }

  \item{term}{
    either a character vector specifying one or more target terms for which nearest neighbours will be found, or a matrix specifying the target vectors directly. A plain vector is interpreted as a single-row matrix.
  }

  \item{n}{
    an integer giving the number of nearest neighbours to be returned for each target term
  }

  \item{M2}{
    an optional dense or sparse matrix (or object of class \code{dsm}). If specified, nearest neighbours are found among the rows (default) or columns (\code{byrow=FALSE}) of \code{M2}, allowing for NN search in a cross-distance setting.
  }
  
  \item{byrow}{
    whether target terms are looked up in rows (default) or columns (\code{byrow=FALSE}) of \code{M}.  NB: Target vectors in the \code{term} argument are always given as row vectors, even if \code{byrow=FALSE}.
  }

  \item{drop}{
    if \code{TRUE}, the return value is simplified to a vector (or distance matrix) if it contains nearest neighbours for exactly one target term (default).  Set \code{drop=FALSE} to ensure that \code{nearest.neighbours} always returns a list.
  }

  \item{skip.missing}{
    if \code{TRUE}, silently ignores target terms not found in the DSM or distance matrix.  By default (\code{skip.missing=FALSE}) an error is raised in this case.
  }

  \item{dist.matrix}{
    if \code{TRUE}, return a full distance matrix between the target term and its nearest neighbours (instead of a vector of neighbours).  Note that a pre-computed distance matrix \code{M} must be symmetric in this case.
  }

  \item{\dots}{
    additional arguments are passed to \code{dist.matrix} if \code{M} is a scored DSM matrix.  See the manpage of \code{\link{dist.matrix}} for details on available parameters and settings.
  }
    
  \item{batchsize}{if \code{term} is a long list of lookup terms, it will automatically be processed in batches.  The number of terms per batch is chosen in such a way that approximately \code{batchsize} intermediate similarity values have to be computed and stored at a time (not used if \code{M} is a pre-computed distance matrix).}
  
  \item{verbose}{if \code{TRUE}, display some progress messages indicating how data are split into batches}
   
}
\details{

  In most cases, the target term itself is automatically excluded from the list of neighbours.  There are two exceptions:
  \enumerate{
    \item{The target term is given as a vector rather than by name.}
    \item{Nearest neighbours are determined in a cross-distance setting.  This is the case if (i) \code{M2} is specified or (ii) \code{M} is a pre-computed distance matrix and not marked to be symmetric.}
  }

  With \code{dist.matrix=TRUE}, the returned distance matrix always includes the target term.
  
  \code{M} can also be a pre-computed distance or similarity matrix from an external source, which must be marked with \code{\link{as.distmat}}.  If \code{M} is a sparse similarity matrix, only non-zero cells will be considered when looking for the nearest neighbours.  Keep in mind that \code{dist.matrix=TRUE} is only valid if \code{M} is a symmetric matrix and marked as such.
  
}
\value{

  A list with one entry for each target \code{term} found in \code{M}, giving
  \itemize{
    \item{\code{dist.matrix=FALSE} (default): the nearest neighbours as a numeric vector of distances or similarities labelled with the corresponding terms and ordered by distance}
    \item{\code{dist.matrix=TRUE}: a full distance or similarity matrix for the target term and its nearest neighbours (as an object of class \code{dist.matrix}). An additional attribute \code{selected} contains a logical vector indicating the position of the target term in the matrix.}
  }

  If \code{drop=TRUE}, a list containing only a single target term will be simplified to a plain vector or distance matrix.

}
\seealso{

  \code{\link{dist.matrix}} for more information on available distance metrics and similarity measures
  
}
\examples{

nearest.neighbours(DSM_Vectors, c("apple_N", "walk_V"), n=10)

nearest.neighbours(DSM_Vectors, "apple_N", n=10, method="maximum")

as.dist(nearest.neighbours(DSM_Vectors, "apple_N", n=10, dist.matrix=TRUE))

}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
