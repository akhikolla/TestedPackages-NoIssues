\name{DSM_Vectors}
\alias{DSM_Vectors}
\docType{data}
\title{
  Pre-Compiled DSM Vectors for Selected Words (wordspace)
}
\description{

  A matrix of 50-dimensional pre-compiled DSM vectors for selected English content words, covering most of the words needed for several basic evaluation tasks included in the package.
  Targets are given as disambiguated lemmas in the form \code{<headword>_<pos>}, e.g. \code{walk_V} and \code{walk_N}.

}
\usage{
  
DSM_Vectors

}
\format{
  A numeric matrix with 1667 rows and 50 columns.
  
  Row labels are disambiguated lemmas of the form \code{<headword>_<pos>}, where the part-of-speech code is one of
  \code{N} (noun), \code{V} (verb), \code{J} (adjective) or \code{R} (adverb).
  
  Attribute \code{"sigma"} contains singular values that can be used for post-hoc power scaling of the latent dimensions (see \code{\link{dsm.projection}}).
}
\details{

  The vocabulary of this DSM covers several basic evaluation tasks, including \code{\link{RG65}}, \code{\link{WordSim353}} and \code{\link{ESSLLI08_Nouns}}, as well as the target nouns \emph{bank} and \emph{vessel} from \code{\link{SemCorWSD}}.  In addition, 40 nearest neighbours each of the words \code{white_J}, \code{apple_N}, \code{kindness_N} and \code{walk_V} are included.

  Co-occurrence frequency data were extracted from a collection of Web corpora with a total size of ca. 9 billion words, using a L4/R4 surface window and 30,000 lexical words as feature terms.  They were scored with sparse simple log-likelihood with an additional log transformation, normalized to Euclidean unit length, and projected into 1000 latent dimensions using randomized SVD (see \code{\link{rsvd}}.  For size reasons, the vectors have been compressed into 50 latent dimensions and renormalized.
  
}
\examples{

nearest.neighbours(DSM_Vectors, "walk_V", 25)

eval.similarity.correlation(RG65, DSM_Vectors) # fairly good

# post-hoc power scaling: whitening (correspond to power=0 in dsm.projection)
sigma <- attr(DSM_Vectors, "sigma")
M <- scaleMargins(DSM_Vectors, cols=1 / sigma)
eval.similarity.correlation(RG65, M) # very good

}
\keyword{datasets}
