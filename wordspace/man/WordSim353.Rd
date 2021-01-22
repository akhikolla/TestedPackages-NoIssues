\name{WordSim353}
\alias{WordSim353}
\docType{data}
\title{
  Similarity Ratings for 351 Noun Pairs (wordspace)
}
\description{

  A database of human similarity ratings for 351 English noun pairs,
  collected by Finkelstein et al. (2002) and annotated with semantic relations
  (similarity vs. relatedness) by Agirre et al. (2009).
  
}
\usage{

WordSim353

}
\format{

  A data frame with 351 rows and the following 6 columns:

  \describe{
    \item{\code{word1}}{first noun (character)}
    \item{\code{word2}}{second noun (character)}
    \item{\code{score}}{average similarity rating by human judges on scale from 0 to 10 (numeric)}
    \item{\code{relation}}{semantic relation between first and second word (factor, see Details below)}
    \item{\code{similarity}}{whether word pair belongs to the \emph{similarity} subset (logical)}
    \item{\code{relatedness}}{whether word pair belongs to the \emph{relatedness} subset (logical)}
  }

  The nouns are given as disambiguated lemmas in the form \code{<headword>_N}.

}
\details{

  The data set is known as \code{WordSim353} because it originally consisted of 353 noun pairs.
  One duplicate entry (\emph{money}--\emph{cash}) as well as the trivial combination 
  \emph{tiger}--\emph{tiger} (which may have been included as a control item)
  have been omitted in the present version, however.

  The following semantic relations are distinguished in the \code{relation} variable:
  \code{synonym}, \code{antonym}, \code{hypernym}, \code{hyponym}, \code{co-hyponym},
  \code{holonym}, \code{meronym} and \code{other} (topically related or completely unrelated).

  Note that the \emph{similarity} and \emph{relatedness} subsets are not disjoint, because they
  share 103 unrelated noun pairs (semantic relation \code{other} and score below 5.0).
}
\source{

  Similarity ratings (Finkelstein \emph{et al.} 2002): \url{http://www.cs.technion.ac.il/~gabr/resources/data/wordsim353/}

  Semantic relations (Agirre \emph{et al.} 2009): \url{http://alfonseca.org/eng/research/wordsim353.html}
  
}
\references{

Agirre, Eneko, Alfonseca, Enrique, Hall, Keith, Kravalova, Jana, Pasca, Marius, and Soroa, Aitor (2009).
A study on similarity and relatedness using distributional and WordNet-based approaches.
In \emph{Proceedings of Human Language Technologies: The 2009 Annual Conference of the North American Chapter of the Association for Computational Linguistics (NAACL-HLT 2009)}, pages 19--27, Boulder, Colorado.

Finkelstein, Lev, Gabrilovich, Evgeniy, Matias, Yossi, Rivlin, Ehud, Solan, Zach, Wolfman, Gadi, and Ruppin, Eytan (2002).
Placing search in context: The concept revisited.
\emph{ACM Transactions on Information Systems}, \bold{20}(1), 116--131.

}
\examples{

head(WordSim353, 20)

table(WordSim353$relation) # semantic relations

# split into "similarity" and "relatedness" subsets
xtabs(~ similarity + relatedness, data=WordSim353) 

}
\keyword{datasets}
