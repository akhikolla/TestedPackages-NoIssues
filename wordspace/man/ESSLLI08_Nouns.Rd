\name{ESSLLI08_Nouns}
\alias{ESSLLI08_Nouns}
\docType{data}
\title{
  Noun Clustering Task from ESSLLI 2008 (wordspace)
}
\description{

  A set of 44 nouns denoting basic-level concepts from 6 semantic classes,
  used as a gold standard in the ESSLLI 2008 shared task on noun clustering.

}
\usage{
  
ESSLLI08_Nouns
  
}
\format{

  A data frame with 44 rows and the following 5 columns:  
  \describe{
    \item{\code{word}}{a character vector specifying the 44 nouns in CWB/Penn format (see \code{\link{convert.lemma}})}
    \item{\code{class}}{a factor vector specifying the semantic class of each noun (bird, fruitTree, green, groundAnimal, tool, vehicle)}
    \item{\code{class2}}{a factor vector specifying a coarser 3-class categorization (animal, vegetable, artifact)}
    \item{\code{class3}}{a factor vector specifying a coarser 2-class categorization (natural, artifact)}
    \item{\code{freq.bnc}}{a numeric vector specifying the frequency of each noun in the British National Corpus}
  }

}
\source{

  \url{http://wordspace.collocations.de/doku.php/data:esslli2008:concrete_nouns_categorization}
  
}
\examples{

print(ESSLLI08_Nouns)

}
\keyword{datasets}
