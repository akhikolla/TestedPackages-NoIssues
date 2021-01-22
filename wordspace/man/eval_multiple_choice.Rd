\name{eval.multiple.choice}
\alias{eval.multiple.choice}
\title{

  Evaluate DSM on Multiple Choice Task (wordspace)

}
\description{

  Evaluates DSM on a multiple choice task by selecting the answer option closest to the target term in distributional space.
  A typical example is the TOEFL Synonym Task (Landauer & Dumais 1997).
  
}
\usage{

eval.multiple.choice(task, M, dist.fnc = pair.distances, \dots,
                     details = FALSE, format = NA, taskname = NA,
                     target.name = "target", correct.name = "correct",
                     distractor.name = "^distract") 

}
\arguments{

  \item{task}{a data frame listing the target word, the correct answer, and one or more additional choices (distractors) for each test item}

  \item{M}{a scored DSM matrix, passed to \code{dist.fnc}}
  
  \item{dist.fnc}{a callback function used to compute distances between term pairs (or similarity scores, which must be marked with an attribute \code{similarity=TRUE}). See \dQuote{Details} below for further information.}

  \item{\dots}{any further arguments are passed to \code{dist.fnc} and can be used e.g. to select a distance measure}

  \item{details}{if \code{TRUE}, a detailed report with information on each task item is returned (see \dQuote{Value} below for details)}
  
  \item{format}{if the task definition specifies POS-disambiguated lemmas in CWB/Penn format, they can automatically be transformed into some other notation conventions; see \code{\link{convert.lemma}} for details}
  
  \item{taskname}{optional row label for the short report (\code{details=FALSE})}

  \item{target.name}{the name of the column of \code{task} containing the target word}

  \item{correct.name}{the name of the column of \code{task} containing the correct choice}
  
  \item{distractor.name}{a regular expression matching columns of \code{task} containing the distractors. The regular expression is matched with \code{perl=TRUE}.}

}
\details{
  
  The callback function \code{dist.fnc} will be invoked with character vectors containing the components of the term pairs as first and second argument,
  the DSM matrix \code{M} as third argument, plus any additional arguments (\code{\dots}) passed to \code{eval.multiple.choice}.
  The return value must be a numeric vector of appropriate length.  If one of the terms in a pair is not represented in the DSM,
  the corresponding distance value should be set to \code{Inf} (or \code{-Inf} in the case of similarity scores).
  In most cases, the default callback \code{\link{pair.distances}} is sufficient if used with suitable parameter settings.

  For each task item, distances between the target word and the possible choices are computed.  Then all choices
  are ranked according to their distances; in the case of a tie, the \emph{higher} rank is assigned to both words.
  A task item counts as a TP (\emph{true positive}, i.e. a successful answer by the DSM) if the correct choice is
  ranked in first place.  Note that if it is tied with another choice, both will be assigned rank 2, so the item does not count as a TP.

  If either the target word is missing from the DSM or none of the choices is found in the DSM, the result for this
  item is set to \code{NA}, which counts as a FP (\emph{false positive}) in the accuracy computation.
  
  With the default \code{dist.fnc} callback, additional arguments \code{method} and \code{p} can be used to select 
  a distance measure (see \code{\link{dist.matrix}} for details).  It is pointless to specify \code{rank="fwd"}, as
  the neighbour ranks produce exactly the same candidate ranking as the distance values.

}
\value{

  The default short report (\code{details=FALSE}) is a data frame with a single row and the columns
  \code{accuracy} (percentage correct), \code{TP} (number of correct answers), \code{FP} (number of wrong answers)
  and \code{missing} (number of test items for which the distance between target and correct choice
  was not found in the DSM).

  The detailed report (\code{details=TRUE}) is a data frame with one row for each task item and the following columns:

  \item{target}{the target word (character)}
  \item{correct}{whether model's choice is correct (logical or \code{NA})}
  \item{best.choice}{best choice according to the DSM (character)}
  \item{best.dist}{distance of best choice from target (numeric)}
  \item{correct.choice}{correct answer (numeric)}
  \item{correct.rank}{rank of correct answer among choices (integer)}
  \item{correct.dist}{distance of correct answer from target (numeric)}

}
\references{

  Landauer, Thomas K. and Dumais, Susan T. (1997).
  A solution to Plato's problem: The latent semantic analysis theory of acquisition, induction and representation of knowledge.
  \emph{Psychological Review}, \bold{104}(2), 211--240.

}
\seealso{

  Suitable gold standard data sets in this package: \bold{TODO}
   
  Support functions: \code{\link{pair.distances}}, \code{\link{convert.lemma}}

}
\examples{
## TODO
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
