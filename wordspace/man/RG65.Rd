\name{RG65}
\alias{RG65}
\docType{data}
\title{
  Similarity Ratings for 65 Noun Pairs (wordspace)
}
\description{

  A database of human similarity ratings for 65 English noun pairs,
  collected by Rubenstein & Goodenough (1965).

}
\usage{

RG65

}
\format{

  A data frame with 65 rows and the following 3 columns:

  \describe{
    \item{\code{word1}}{first noun (character)}
    \item{\code{word2}}{second noun (character)}
    \item{\code{score}}{average similarity rating by human judges on scale from 0 to 4 (numeric)}
  }

  The nouns are given as disambiguated lemmas in the form \code{<headword>_N}.

}
\details{

  The word pairs are sorted by increasing similarity score, as in the original paper.

}
\source{

  Rubenstein, Herbert and Goodenough, John B. (1965).
  Contextual correlates of synonymy.
  \emph{Communications of the ACM}, \bold{8}(10), 627--633.

}
%% \references{
%%  ~~ perhaps references to evaluation studies using RG65 ~~
%% }
\examples{

head(RG65, 10) # least similar pairs
tail(RG65, 10) # most similar pairs

}
\keyword{datasets}
