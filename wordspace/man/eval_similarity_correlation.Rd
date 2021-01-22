\name{eval.similarity.correlation}
\alias{eval.similarity.correlation}
\title{

  Evaluate DSM on Correlation with Similarity Ratings (wordspace)

}
\description{

  Performs evaluation by comparing the distances (or similarities) computed by a DSM with (typically human) word similarity ratings.
  Well-know examples are the noun pair ratings collected by Rubenstein & Goodenough (1965; \code{\link{RG65}}) and Finkelstein et al. (2002;  \code{\link{WordSim353}}).
  
  The quality of the DSM predictions is measured by Spearman rank correlation \eqn{rho}.
  
}
\usage{
eval.similarity.correlation(task, M, dist.fnc=pair.distances,
                            details=FALSE, format=NA, taskname=NA,
                            word1.name="word1", word2.name="word2", score.name="score",
                            \dots)
}
\arguments{

  \item{task}{a data frame containing word pairs (usually in columns \code{word1} and \code{word2}) with similarity ratings (usually in column \code{score}); any other columns will be ignored}

  \item{M}{a scored DSM matrix, passed to \code{dist.fnc}}
  
  \item{dist.fnc}{a callback function used to compute distances or similarities between word pairs.
    It will be invoked with character vectors containing the components of the word pairs as first and second argument,
    the DSM matrix \code{M} as third argument, plus any additional arguments (\code{\dots}) passed to \code{eval.similarity.correlation}.
    The return value must be a numeric vector of appropriate length.  If one of the words in a pair is not represented in the DSM,
    the corresponding distance value should be set to \code{Inf} (or \code{-Inf} in the case of similarities).
  }

  \item{details}{if \code{TRUE}, a detailed report with information on each task item is returned (see Value below for details)}
  
  \item{format}{if the task definition specifies POS-disambiguated lemmas in CWB/Penn format, they can automatically be transformed into some other notation conventions; see \code{\link{convert.lemma}} for details}
  
  \item{taskname}{optional row label for the short report (\code{details=FALSE})}
  
  \item{\dots}{any further arguments are passed to \code{dist.fnc} and can be used e.g. to select a distance measure}

  \item{word1.name}{the name of the column of \code{task} containing the first word of each pair}

  \item{word2.name}{the name of the column of \code{task} containing the second word of each pair}

  \item{score.name}{the name of the column of \code{task} containing the corresponding similarity ratings}

}
\details{
  
  DSM distances are computed for all word pairs and compared with similarity ratings from the gold standard.
  As an evaluation criterion, Spearman rank correlation between the DSM and gold standard scores is computed.
  The function also reports a confidence interval for Pearson correlation, which might require suitable transformation
  to ensure a near-linear relationship in order to be meaningful.
  
  \bold{NB:} Since the correlation between similarity ratings and DSM distances will usually be negative, the evaluation
  report omits minus signs on the correlation coefficients.

  With the default \code{dist.fnc}, the distance values can optionally be transformed through an arbitrary function specified in the \code{transform} argument (see \code{\link{pair.distances}} for details).
  Examples include \code{transform=log} (esp. for neighbour rank as a distance measure) 
  and \code{transform=function (x) 1/(1+x)} (in order to transform distances into similarities).
  Note that Spearman rank correlation is not affected by any monotonic transformation, so the main evaluation results
  will remain unchanged.

  If one or both words of a pair are not found in the DSM, the distance is set to a fixed value 10\% above the
  maximum of all other DSM distances, or 10\% below the minimum in the case of similarity values.
  This is done in order to avoid numerical and visualization problems with \code{Inf} values;
  the particular value used does not affect the rank correlation coefficient.
    
  With the default \code{dist.fnc} callback, additional arguments \code{method} and \code{p} can be used to select 
  a distance measure (see \code{\link{dist.matrix}} for details); \code{rank=TRUE} can be specified in order to 
  use neighbour rank as a measure of semantic distance.

}
\value{

  The default short report (\code{details=FALSE}) is a data frame with a single row and the following columns:

  \item{rho}{(absolute value of) Spearman rank correlation coefficient \eqn{rho}}
  \item{p.value}{p-value indicating evidence for a significant correlation}
  \item{missing}{number of pairs not included in the DSM}
  \item{r}{(absolute value of) Pearson correlation coefficient \eqn{r}}
  \item{r.lower}{lower bound of confidence interval for Pearson correlation}
  \item{r.upper}{upper bound of confidence interval for Pearson correlation}

  The detailed report (\code{details=TRUE}) is a copy of the original task data with two additional columns:

  \item{distance}{distance calculated by the DSM for each word pair, possibly transformed (numeric)}
  \item{missing}{whether word pair is missing from the DSM (logical)}

  In addition, the short report is appended to the data frame as an attribute \code{"eval.result"}, 
  and the optional \code{taskname} value as attribute \code{"taskname"}.  The data frame is marked as an
  object of class \code{eval.similarity.correlation}, for which suitable \code{\link[=print.eval.similarity.correlation]{print}}
  and \code{\link[=plot.eval.similarity.correlation]{plot}} methods are defined.

}
\references{

Finkelstein, Lev, Gabrilovich, Evgeniy, Matias, Yossi, Rivlin, Ehud, Solan, Zach, Wolfman, Gadi, and Ruppin, Eytan (2002).
Placing search in context: The concept revisited.
\emph{ACM Transactions on Information Systems}, \bold{20}(1), 116--131.

Rubenstein, Herbert and Goodenough, John B. (1965).
Contextual correlates of synonymy.
\emph{Communications of the ACM}, \bold{8}(10), 627--633.

}
\seealso{
  Suitable gold standard data sets in this package: \code{\link{RG65}}, \code{\link{WordSim353}}
  
  Support functions: \code{\link{pair.distances}}, \code{\link{convert.lemma}}
  
  Plotting and printing evaluation results: \code{\link{plot.eval.similarity.correlation}}, \code{\link{print.eval.similarity.correlation}}
}
\examples{

eval.similarity.correlation(RG65, DSM_Vectors)

\dontrun{
plot(eval.similarity.correlation(RG65, DSM_Vectors, details=TRUE))
}
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
