\name{SemCorWSD}
\alias{SemCorWSD}
\encoding{UTF-8}
\docType{data}
\title{
  SemCor Word Sense Disambiguation Task (wordspace)
}
\description{

  A collection of sentences containing ambiguous words manually labelled with WordNet senses.
  The data were obtained from the SemCor corpus version 3.0.

}
\usage{

SemCorWSD

}
\format{

  A data frame with 647 rows and the following 8 columns (all of type character):

  \describe{
    \item{\code{id}}{Unique item ID}
    \item{\code{target}}{The target word (lemmatized)}
    \item{\code{pos}}{Word class of the target word (\code{n}, \code{v} or \code{a})}
    \item{\code{sense}}{Sense of the target word in this sentence (given as a WordNet lemma)}
    \item{\code{gloss}}{WordNet definition of this sense}
    \item{\code{sentence}}{The sentence containing the ambiguous word}
    \item{\code{hw}}{Lemmatized form of the sentence (\dQuote{headwords}); punctuation marks are excluded and all remaining words are case-folded}
    \item{\code{lemma}}{Lemmatized and POS-disambiguated form in CWB/Penn format, e.g. \code{move_N} for the headword \emph{move} used as a noun}
  }

}
\details{

  Target words and senses had to meet the following criteria in order to be included in the data set:
  \itemize{
    \item{sense occurs \eqn{f \ge 5} times in SemCor 3.0}
    \item{sense accounts for at least 10\% of all occurrences of the target}
    \item{at least two senses of target remain after previous two filters}
  }

  \code{SemCorWSD} contains sentence contexts for the following target words:
  \itemize{
    \item{ambiguous nouns from Schütze (1998): \emph{interest}, \emph{plant}, \emph{space}, \emph{vessel}}
    \item{misc. ambiguous nouns: \emph{bank}}
    \item{misc. ambiguous verbs: \emph{find}, \emph{grasp}, \emph{open}, \emph{run}}
  }
  
}
\source{

  \bold{TODO} (SemCor reference, NLTK extraction)

}
\seealso{

  \code{\link{context.vectors}}
  
}
\references{

  Schütze, Hinrich (1998). Automatic word sense discrimination. \emph{Computational Linguistics}, \bold{24}(1), 97--123.

}
\examples{

with(SemCorWSD, table(sense, target))

# all word senses with brief definitions ("glosses")
with(SemCorWSD, sort(unique(paste0(target, " ", sense, ": ", gloss))))

}
\keyword{datasets}
