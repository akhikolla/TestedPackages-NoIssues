\name{pair.distances}
\alias{pair.distances}
\title{

  Semantic Distances Between Word Pairs (wordspace)

}
\description{

  Compute semantic distances (or similarities) between pairs of target terms based on a scored DSM matrix \code{M},
  according to any of the distance measures supported by \code{\link{dist.matrix}}.
  If one of the terms in a pair is not represented in the DSM, the distance is set to \code{Inf}
  (or to \code{-Inf} in the case of a similarity measure).

}
\usage{

pair.distances(w1, w2, M, \dots, transform = NULL, 
               rank = c("none", "fwd", "bwd", "avg"),
               avg.method = c("arithmetic", "geometric", "harmonic"),
               batchsize = 10e6, verbose = FALSE)

}
\arguments{

  \item{w1}{a character vector specifying the first term of each pair}

  \item{w2}{a character vector of the same length as \code{w1}, specifying the second term of each pair}
  
  \item{M}{a sparse or dense DSM matrix, suitable for passing to \code{\link{dist.matrix}}, or an object of class \code{dsm}. Alternatively, \code{M} can be a pre-computed distance or similarity matrix returned by \code{\link{dist.matrix}} or marked as such with \code{\link{as.distmat}}.}
  
  \item{\dots}{further arguments are passed to \code{dist.matrix} and determine the distance or similarity measure to be used (see \code{\link{dist.matrix}} for details)}

  \item{rank}{whether to return the distance between the two terms (\code{"none"}) or the neighbour rank (see \dQuote{Details} below)}

  \item{transform}{an optional transformation function applied to the distance, similarity or rank values (e.g. \code{transform=log10} for logarithmic ranks). This option is provided as a convenience for evaluation code that calls \code{pair.distances} with user-specified arguments.}

  \item{avg.method}{with \code{rank="avg"}, whether to compute the arithmetic, geometric or harmonic mean of forward and backward rank}
  
  \item{batchsize}{maximum number of similarity values to compute per batch. This parameter has an essential influence on efficiency and memory use of the algorithm and has to be tuned carefully for optimal performance.}
  
  \item{verbose}{if \code{TRUE}, display some progress messages indicating how data are split into batches}

}
\details{

  The \code{rank} argument controls whether semantic distance is measured directly by geometric distance (\code{none}),
  by forward neighbour rank (\code{fwd}), by backward neighbour rank (\code{bwd}), or by the average of forward and backward rank (\code{avg}).
  Forward neighbour rank is the rank of \code{w2} among the nearest neighbours of \code{w1}.
  Backward neighbour rank is the rank of \code{w1} among the nearest neighbours of \code{w2}.
  The average can be computed as an arithmetic, geometric or harmonic mean, depending on \code{avg.method}.
  
  Note that a transformation function is applied \emph{after} averaging.
  In order to compute the arithmetic mean of log ranks, set \code{transform=log10}, \code{rank="avg"} and \code{avg.method="geometric"}.

  Neighbour ranks assume that each target term is its own nearest neighbour and adjust ranks to account for this (i.e. \code{w1 == w2} should return a rank of 0).
  If \code{M} is a pre-computed distance matrix, the adjustment is only applied if it is also marked as symmetric (because otherwise \code{w1} might not appear in the list of neighbours at all).  This might lead to unexpected results once asymmetric measures are implemented in \code{dist.matrix}.
  
  For a sparse pre-computed similarity matrix \code{M}, only non-zero cells are considered as neighbours and all other ranks are set to \code{Inf}.  This is consistent with the behaviour of \code{\link{nearest.neighbours}}.

  \code{pair.distances} is used as a default callback in several evaluation functions, which rely on the attribute \code{similarity} to distinguish between distance measures and similarity scores.  For this reason, transformation functions should always be \bold{isotonic} (order-preserving) so as not to mislead the evaluation procedure.

}
\value{
  
  If \code{rank="none"} (the default), a numeric vector of the same length as \code{w1} and \code{w2}
  specifying the distances or similarities between the term pairs, according to the metric selected with the extra arguments (\code{\dots}).
  
  Otherwise, an integer or numeric vector of the same length as \code{w1} and \code{w2} specifying
  forward, backward or average neighbour rank for the two terms.
  
  In either case, a distance or rank of \code{Inf} (or a similarity of \code{-Inf}) is returned for any term pair not represented in the DSM.
  Attribute \code{similarity} is set to \code{TRUE} if the returned values are similarity scores rather than distances.

}
\seealso{
  \code{\link{dist.matrix}}, \code{\link{eval.similarity.correlation}}, \code{\link{eval.multiple.choice}}, \code{\link{nearest.neighbours}}
}
\examples{

transform(RG65, angle=pair.distances(word1, word2, DSM_Vectors))

}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
