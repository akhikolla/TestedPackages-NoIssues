\name{context.vectors}
\alias{context.vectors}
\encoding{utf8}
\title{
  Compute Bag-of-Words Context Vectors (wordspace)
}
\description{

  Compute bag-of-words context vectors as proposed by Schütze (1998) for automatic word sense disambiguation and induction.  Each context vector is the centroid of the DSM vectors of all terms occurring in the context.

}
\usage{

context.vectors(M, contexts, split = "\\\\s+",
                drop.missing = TRUE, row.names=NULL)

}
\arguments{

  \item{M}{
    numeric matrix of row vectors for the terms specified by \code{rownames(M)}, or an object of class \code{dsm}
  }

  \item{contexts}{
    the contexts for which bag-of-words representations are to be computed.  Must be a character vector, a list of character vectors, or a list of labelled numeric vectors (see Details below).
  }

  \item{split}{
    Perl regular expression determining how contexts given as a character vector are split into terms.  The default behaviour is to split on whitespace.
  }

  \item{drop.missing}{
    if \code{TRUE} (default), contexts that do not contain any known terms are silently dropped; otherwise the corresponding context vectors will be all zeroes.
  }

  \item{row.names}{
    a character vector of the same length as \code{contexts}, specifying row names for the resulting matrix of centroid vectors
  }
}
\details{

  The \code{contexts} argument can be specified in several different ways:
  \itemize{
    \item A character vector: each element represents a context given as a string, which will be split on the Perl regular expression \code{split} and then looked up in \code{M}. Repetitions are allowed and will be weighted accordingly in the centroid.
    \item A list of character vectors: each item represents a pre-tokenized context given as a sequence of terms to be looked up in \code{M}. Repetitions are allowed and will be weighted accordingly in the centroid.
    \item A list of labelled numeric vectors: each item represents a bag-of-words representation of a context, where labels are terms to be looked up in \code{M} and the corresponding values their frequency counts or (possibly non-integer) weights.
    \item \emph{(deprecated)} A logical vector corresponding to the rows of \code{M}, which will be used directly as an index into \code{M}.
    \item \emph{(deprecated)} An unlabelled integer vector, which will be used as an index into the rows of \code{M}.
  }
  
  For each context, terms not found in the matrix \code{M} are silently computed.  Then a context vector is computed as the centroid of the remaining term vectors.  If the context contains multiple occurrences of the same term, its vector will be weighted accordingly.  If the context is specified as a bag-of-words representations, the terms are weighted according to the corresponding numerical values.
  
  Neither word order nor any other structural properties of the contexts are taken into account.
}
\value{

  A numeric matrix with the same number of columns as \code{M} and one row for each context (excluding contexts without known terms if \code{drop.missing=TRUE}).  If the vector \code{contexts} has names or \code{row.names} is specified, the matrix rows will be labelled accordingly.  Otherwise the row labels correspond to the indices of the respective entries in \code{contexts}, so matrix rows can always be identified unambiguously if \code{drop.missing=TRUE}.

  If \code{drop.missing=FALSE}, a context without any known terms (including an empty context) is represented by an all-zero vector.
  
}
\references{

  Schütze, Hinrich (1998). Automatic word sense discrimination. \emph{Computational Linguistics}, \bold{24}(1), 97--123.  

}
\seealso{

  \code{\link{SemCorWSD}}
  
}
\examples{
# different ways of specifying contexts
M <- DSM_TermTermMatrix
context.vectors(M, c("dog cat cat", "cause effect")) # contexts as strings
context.vectors(M, list(c("dog", "cat", "cat"), c("cause", "effect"))) # pre-tokenized
context.vectors(M, list(c(dog=1, cat=2), c(cause=1, effect=1))) # bag of words

# illustration of WSD algorithm: 6 sentences each for two senses of "vessel"
VesselWSD <- subset(SemCorWSD, target == "vessel")
with(VesselWSD, cat(paste0(sense, ": ", sentence, "\n")))

# provide sense labels in case some contexts are dropped b/c of too many missing words
Centroids <- with(VesselWSD, context.vectors(DSM_Vectors, lemma, row.names=sense))
Centroids[, 1:5]

(res <- kmeans(Centroids, 2)$cluster) # flat clustering with k-means
table(rownames(Centroids), res)       # ... works perfectly

\dontrun{
plot(hclust(dist.matrix(Centroids, as.dist=TRUE)))
}
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
