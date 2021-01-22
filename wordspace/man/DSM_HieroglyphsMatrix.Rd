\name{DSM_HieroglyphsMatrix}
\alias{DSM_HieroglyphsMatrix}
\docType{data}
\title{
  A Small Co-occurrence Matrix (wordspace)
}
\description{

  A small co-occurrence matrix of verb-object combinations from the British National Corpus (BNC)
  Verbs correspond to columns of the matrix and their object nouns to rows.
  This matrix is shown as the "hieroglyphs" example in the DSM turorial.

}
\usage{

DSM_HieroglyphsMatrix

}
\format{
  
  A numeric matrix with 7 rows and 6 columns.

  Rows represent the target nouns \emph{knife}, \emph{cat}, \emph{dog}, \emph{boat}, \emph{cup}, \emph{pig} and \emph{banana}.
  Columns represent the feature verbs \emph{get}, \emph{see}, \emph{use}, \emph{hear}, \emph{eat} and \emph{kill}.
  
}
%% \details{
%% }
%% \source{
%% }
\examples{

print(DSM_HieroglyphsMatrix)

## cosine similarities between rows of the matrix
round(dist.matrix(DSM_HieroglyphsMatrix, convert=FALSE), 3)

}
\keyword{datasets}
