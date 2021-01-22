\name{DSM_GoodsMatrix}
\alias{DSM_GoodsMatrix}
\docType{data}
\title{
  A Scored Co-occurrence Matrix of Nouns Denoting Goods (wordspace)
}
\description{

  A pre-scored verb-object co-occurrence matrix for 240 target nouns denoting goods
  and the 3 feature verbs \emph{own}, \emph{buy} and \emph{sell}.
  This matrix is useful for illustrating the application and purpose of dimensionality reduction techniques.

}
\usage{

DSM_GoodsMatrix

}
\format{
  
  A numeric matrix with 240 rows corresponding to target nouns denoting goods and 4 columns, corresponding to
  \describe{
    \item{\code{own}, \code{buy}, \code{sell}:}{
      association scores for co-occurrences of the nouns with the verbs \emph{own}, \emph{buy} and \emph{sell}
    }
    \item{\code{fringe}:}{
      an indicator of how close each point is to the \dQuote{fringe} of the data set (ranging from 0 to 1)
    }
  }

}
\details{

  Co-occurrence data are based on verb-object dependency relations in the British National Corpus, obtained from \code{\link{DSM_VerbNounTriples_BNC}}.  Only nouns that co-occur with all three verbs are included in the data set.

  The co-occurrence matrix is weighted with \emph{non-sparse} log-likelihood (\code{simple-ll}) and an additional logarithmic transformation (\code{log}).  Row vectors are \emph{not} normalized.

  The \emph{fringeness score} in column \code{fringe} indicates how close a data point is to the fringe of the data set.  Values are distance quantiles based on PCA-whitened Manhattan distance from the centroid.  For example, \code{fringe >= .8} characterizes 20\% of points that are closest to the fringe.  Fringeness is mainly used to select points to be labelled in plots or to take stratified samples from the data set.
  
}
\examples{

DSM_GoodsMatrix[c("time", "goods", "service"), ]

}
\keyword{datasets}
