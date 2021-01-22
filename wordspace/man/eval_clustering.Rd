\name{eval.clustering}
\alias{eval.clustering}
\title{

  Evaluate DSM on Clustering Task (wordspace)

}
\description{

  Performs evaluation on a word clustering task by comparing a flat clustering
  solution based on semantic distances with a gold classification. 
  
}
\usage{

eval.clustering(task, M, dist.fnc = pair.distances, \dots,
                details = FALSE, format = NA, taskname = NA,
                scale.entropy = FALSE, n.clusters = NA,
                word.name = "word", class.name = "class")

}
\arguments{

  \item{task}{a data frame listing words and their classes, usually in columns named \code{word} and \code{class}}

  \item{M}{a scored DSM matrix, passed to \code{dist.fnc}}

  \item{dist.fnc}{a callback function used to compute distances between word pairs.
    It will be invoked with character vectors containing the components of the word pairs as first and second argument,
    the DSM matrix \code{M} as third argument, plus any additional arguments (\code{\dots}) passed to \code{eval.multiple.choice}.
    The return value must be a numeric vector of appropriate length.  If one of the words in a pair is not represented in the DSM,
    the corresponding distance value should be set to \code{Inf}.
  }

  \item{\dots}{any further arguments are passed to \code{dist.fnc} and can be used e.g. to select a distance measure}
    
  \item{details}{if \code{TRUE}, a detailed report with information on each task item is returned (see \dQuote{Value} below for details)}
  
  \item{format}{if the task definition specifies POS-disambiguated lemmas in CWB/Penn format, they can automatically be transformed into some other notation conventions; see \code{\link{convert.lemma}} for details}
  
  \item{taskname}{optional row label for the short report (\code{details=FALSE})}

  \item{scale.entropy}{whether to scale cluster entropy values to the range \eqn{[0, 1]}}

  \item{n.clusters}{number of clusters. The (very sensible) default is to generate as many clusters as their are classes in the gold standard.}
  
  \item{word.name}{the name of the column of \code{task} containing words}

  \item{class.name}{the name of the column of \code{task} containing gold standard classes}

}
\details{

  The test words are clustered using the \dQuote{partitioning around medoids} (PAM) algorithm (Kaufman \& Rousseeuw 1990, Ch. 2) based on their semantic distances.
  The PAM algorithm is used because it works with arbitrary distance measures (including neihbour rank), produces a stable solution (unlike most iterative algorithms)
  and has shown to be on par with state-of-the-art spherical k-means clustering (CLUTO) in evaluation studies.

  Each cluster is automatically assigned a majority label, i.e. the gold standard class occurring most frequently in the cluster.
  This represents the best possible classification that can be derived from the clustering.

  As evaluation metrics, clustering \bold{purity} (accuracy of the majority classification) and \bold{entropy} are computed.
  The latter is defined as a weighted average over the entropy of the class distribution within each cluster, expressed in bits.
  If \code{scale.entropy=TRUE}, the value is divided by the overall entropy of the class distribution in the gold standard, scaling it to the range \eqn{[0, 1]}.

  NB: The semantic distance measure selected with the extra arguments (\code{\dots}) should be \emph{symmetric}.
  In particular, it is not very sensible to specify \code{rank="fwd"} or \code{rank="bwd"}.

  NB: Similarity measures are not supported by the current clustering algorithm.  Make sure not to call \code{dist.matrix}
  (from \code{dist.fnc}) with \code{convert=FALSE}!
}
\value{

  The default short report (\code{details=FALSE}) is a data frame with a single row and the columns
  \code{purity} (clustering purity as a percentage), \code{entropy} (scaled or unscaled clustering entropy)
  and \code{missing} (number of words not found in the DSM).

  The detailed report (\code{details=TRUE}) is a data frame with one row for each test word and the following columns:

  \item{word}{the test word (character)}
  \item{cluster}{cluster to which the word has been assigned; all unknown words are collected in an additional cluster \code{"n/a"}}
  \item{label}{majority label of this cluster (factor with same levels as \code{gold})}
  \item{gold}{gold standard class of the test word (factor)}
  \item{correct}{whether majority class assignment is correct (logical)}
  \item{missing}{whether word was not found in the DSM (logical)}

}
\seealso{

  Suitable gold standard data sets in this package: \code{\link{ESSLLI08_Nouns}}
  
  Support functions: \code{\link{pair.distances}}, \code{\link{convert.lemma}}
}
\examples{

eval.clustering(ESSLLI08_Nouns, DSM_Vectors, class.name="class2")

}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
