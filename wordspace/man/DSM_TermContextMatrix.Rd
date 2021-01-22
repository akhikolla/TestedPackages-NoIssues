\name{DSM_TermContextMatrix}
\alias{DSM_TermContextMatrix}
\alias{DSM_TermContext}
\docType{data}
\title{
  Example of a Term-Context Co-occurrence Matrix (wordspace)
}
\description{

  This matrix is a typical example of a term-context DSM co-occurrence matrix,
  derived from the English Wikipedia.  It is available as a plain matrix in sparse representation,
  and as DSM object including marginal frequency data.

}
\usage{

DSM_TermContextMatrix

DSM_TermContext

}
\format{

  \code{DSM_TermContextMatrix} is a sparse numeric matrix of class \code{\link[=dgCMatrix-class]{dgCMatrix}} with 7 rows and 7 columns.

  Rows represent the target nouns \emph{cat}, \emph{dog}, \emph{animal}, \emph{time}, \emph{reason}, \emph{cause}, \emph{effect}.

  Columns specify the occurrence frequencies of these nouns in Wikipedia articles on \emph{Felidae}, \emph{Pet}, \emph{Feral}, \emph{Boat}, \emph{Philosohpy}, \emph{Kant} and \emph{Back Pain}.

  \code{DSM_TermContext} is an object of class \code{\link{dsm}} based on the same co-occurrence matrix, but with additional information on marginal frequencies of the target terms and feature contexts.
  
}
%% \details{
%% }
%% \source{
%%
%% - WP500 corpus
%% - term-document frequencies extracted with UCS toolkit
%% - marginal frequencies based on full matrix, not actual article sizes
%% - example DSM is a subset of WP500_TermDoc_Lemma, with minor edits
%%
%% }
\seealso{

  This matrix/DSM describes the same target nouns as the term-term matrix \code{\link{DSM_TermTermMatrix}} and corresponding DSM object \code{\link{DSM_TermTerm}}.

}
\examples{

DSM_TermContextMatrix["time", ] # row vector for target noun "time"

all.equal(DSM_TermContextMatrix, head(DSM_TermContext, Inf))

# M M' = symmetric matrix of co-occurrence frequencies of nouns within articles
tcrossprod(DSM_TermContextMatrix) 

}
\keyword{datasets}
