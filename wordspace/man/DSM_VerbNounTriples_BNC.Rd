\name{DSM_VerbNounTriples_BNC}
\alias{DSM_VerbNounTriples_BNC}
\docType{data}
\title{
  Verb-Noun Co-occurrence Frequencies from British National Corpus (wordspace)
}
\description{

  A table of co-occurrence frequency counts for verb-subject and verb-object pairs in the British National Corpus (BNC).
  Subject and object are represented by the respective head noun.  Both verb and noun entries are lemmatized.
  Separate frequency counts are provided for the written and the spoken part of the BNC.

}
\usage{
  
DSM_VerbNounTriples_BNC

}
\format{
  A data frame with 250117 rows and the following columns:
  \describe{
    \item{\code{noun}:}{noun lemma}
    \item{\code{rel}:}{syntactic relation (\code{subj} or \code{obj})}
    \item{\code{verb}:}{verb lemma}
    \item{\code{f}:}{co-occurrence frequency of noun-rel-verb triple in subcorpus}
    \item{\code{mode}:}{subcorpus (\code{written} for the writte part of the BNC, \code{spoken} for the spoken part of the BNC)}
  }
  
}
\details{
  In order to save disk space, triples that occur less than 5 times in the respective subcorpus have been omitted from the table.
  The data set should therefore not be used for practical applications.
}
\source{
  Syntactic dependencies were extracted from the British National Corpus (Aston & Burnard 1998) using the C&C robust syntactic parser (Curran \emph{et al.} 2007).  Lemmatization and POS tagging are also based on the C&C output.
}
\references{
  Aston, Guy and Burnard, Lou (1998). \emph{The BNC Handbook}. Edinburgh University Press, Edinburgh. See also the BNC homepage at \url{http://www.natcorp.ox.ac.uk/}.

  Curran, James; Clark, Stephen; Bos, Johan (2007). Linguistically motivated large-scale NLP with {C\&C} and Boxer. In \emph{Proceedings of the 45th Annual Meeting of the Association for Computational Linguistics, Posters and Demonstrations Sessions}, pages 33--36, Prague, Czech Republic.
}
\examples{

# compile some typical DSMs for spoken part of BNC
bncS <- subset(DSM_VerbNounTriples_BNC, mode == "spoken")
dim(bncS) # ca. 14k verb-rel-noun triples

# dependency-filtered DSM for nouns, using verbs as features
# (note that multiple entries for same relation are collapsed automatically)
bncS_depfilt <- dsm(
  target=bncS$noun, feature=bncS$verb, score=bncS$f,
  raw.freq=TRUE, verbose=TRUE)

# dependency-structured DSM
bncS_depstruc <- dsm(
  target=bncS$noun, feature=paste(bncS$rel, bncS$verb, sep=":"), score=bncS$f,
  raw.freq=TRUE, verbose=TRUE)

}
\keyword{datasets}
