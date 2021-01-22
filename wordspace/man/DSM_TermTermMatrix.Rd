\name{DSM_TermTermMatrix}
\alias{DSM_TermTermMatrix}
\alias{DSM_TermTerm}
\docType{data}
\title{
  Example of a Term-Term Co-occurrence Matrix (wordspace)
}
\description{

  This matrix is a typical example of a term-term DSM co-occurrence matrix,
  derived from the English Wikipedia.  It is available as a plain matrix in dense representation,
  and as a DSM object including marginal frequency data.

}
\usage{

DSM_TermTermMatrix

DSM_TermTerm

}
\format{

  \code{DSM_TermTermMatrix} is a numeric matrix with 7 rows and 7 columns.

  Rows represent the target nouns \emph{cat}, \emph{dog}, \emph{animal}, \emph{time}, \emph{reason}, \emph{cause}, \emph{effect}.

  Columns specify co-occurrence frequencies of these nouns with the words \emph{breed}, \emph{tail}, \emph{feed}, \emph{kill}, \emph{important}, \emph{explain} and \emph{likely} in articles of the English Wikipedia.  Co-occurring words must appear within a distance of at most two word tokens of each other.

  \code{DSM_TermTerm} is an object of class \code{\link{dsm}} based on the same co-occurrence matrix, but with additional information on marginal frequencies of the target and feature terms.

}
%% \details{
%% }
%% \source{
%%
%% - WP500 corpus
%% - term-term co-occurrences within L2/R2 span (skipping punctuation, limited to sentences) extracted with UCS toolkit
%% - example DSM is a subset of WP500_TermDoc_Lemma (unedited)
%%
%% }
\seealso{

  This matrix/DSM describes the same target terms as the term-context matrix \code{\link{DSM_TermContextMatrix}} and corresponding DSM object \code{\link{DSM_TermContext}}.

}
\examples{

DSM_TermTermMatrix["time", ] # row vector for target noun "time"

all.equal(DSM_TermTermMatrix, head(DSM_TermTerm, Inf))

\dontrun{
plot(hclust(dist.matrix(DSM_TermTermMatrix, as.dist=TRUE)))
}

}
\keyword{datasets}
