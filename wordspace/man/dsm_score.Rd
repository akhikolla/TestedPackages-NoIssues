\name{dsm.score}
\alias{dsm.score}
\encoding{utf8}
\title{
  Weighting, Scaling and Normalisation of Co-occurrence Matrix (wordspace)
}
\description{

  Compute feature scores for a term-document or term-term co-occurrence matrix, using one of several standard association measures.  Scores can optionally be rescaled with an isotonic transformation function and centered or standardized.  In addition, row vectors can be normalized to unit length wrt. a given norm. 

  This function has been optimized for efficiency and low memory overhead.
  
}
\usage{

dsm.score(model, score = "frequency",
          sparse = TRUE, negative.ok = NA,
          transform = c("none", "log", "root", "sigmoid"),
          scale = c("none", "standardize", "center", "scale"),
          normalize = FALSE, method = "euclidean", p = 2, tol = 1e-6,
          matrix.only = FALSE, update.nnzero = FALSE,
          batchsize = 1e6, gc.iter = Inf)

}
\arguments{
  
  \item{model}{a DSM model, i.e. an object of class \code{dsm}}

  \item{score}{the association measure to be used for feature weighting; either a character string naming one of the built-in measures or a user-defined function (see \dQuote{Details} below)}

  \item{sparse}{if \code{TRUE} (the default), compute sparse non-negative association scores (see \dQuote{Details} below).
    Non-sparse association scores are only allowed if \code{negative.ok=TRUE}.
  }

  \item{negative.ok}{whether operations that introduce negative values into the score matrix (non-sparse association scores, standardization of columns, etc.) are allowed.
    The default (\code{negative.ok=NA}) is \code{TRUE} if the co-occurrence matrix \eqn{M} is dense, and \code{FALSE} if it is sparse.  See \dQuote{Details} below for the special value \code{negative.ok="nonzero"}. 
  }
    
  \item{transform}{scale transformation to be applied to association scores (see \dQuote{Details} below)}

  \item{scale}{if not \code{"none"}, standardize columns of the scored matrix by z-transformation (\code{"standardize"}), center them without rescaling (\code{"center"}), or scale to unit RMS without centering (\code{"scale"})}

  \item{normalize}{if \code{TRUE} normalize row vectors of scored matrix to unit length, according to the norm indicated by \code{method} and \code{p}}

  \item{method, p}{norm to be used with \code{normalize=TRUE}.
    See \code{\link{rowNorms}} for admissible values and details on the corresponding norms
  }
  
  \item{tol}{if \code{normalize=TRUE}, row vectors with norm below \code{tol} are explicitly set to all zeroes instead of attempting to normalize them (see \code{\link{normalize.rows}} for more information)}
  
  \item{matrix.only}{whether to return updated DSM model (default) or only the matrix of scores (\code{matrix.only=TRUE})}
  
  \item{update.nnzero}{if \code{TRUE} and a full DSM model is returned, update the counts of nonzero entries in rows and columns according to the matrix of scores (there may be fewer nonzero entries with sparse association scores, or more from dense association scores and/or column scaling)}

  \item{batchsize}{if \code{score} is a user-defined function, the co-occurrence matrix is divided into blocks of approx. \code{batchsize} elements each in order to reduce memory overhead}

  \item{gc.iter}{how often to run the garbage collector when computing user-defined association scores; \code{gc()} is called after every \code{gc.iter} batches in order to reclaim temporary data and keep memory overhead as low as possible. This option should only be specified if memory is very tight, since garbage collector runs can be expensive (e.g. when there are many distinct strings in the workspace). With the default value \code{gc.iter=Inf}, no calls to \code{gc()} will be made; \code{gc.iter=4} seems to give a good trade-off between memory overhead and degraded performance.}
}
\details{
  \subsection{Association measures}{
    Association measures (AM) for feature scoring are defined in the notation of Evert (2008).  The most important symbols are \eqn{O_{11} = O}{O11 = O} for the observed co-occurrence frequency, \eqn{E_{11} = E}{E11 = E} for the co-occurrence frequency expected under a null hypothesis of independence, \eqn{R_1}{R1} for the marginal frequency of the target term, \eqn{C_1}{C1} for the marginal frequency of the feature term or context, and \eqn{N} for the sample size of the underlying corpus.  Evert (2008) explains in detail how these values are computed for different types of co-occurrence; practical examples can be found in the distributional semantics tutorial at \url{http://wordspace.collocations.de/}.

    Several commonly used AMs are implemented in optimized C++ code for efficiency and minimal memory overhead.  They are selected by name, which is passed as a character string in the \code{score} argument.  See below for a list of built-in measures and their full equations.
    
    Other AMs can be applied by passing a user-defined function in the \code{score} argument. See \dQuote{User-defined association measures} at the end of this section for details.
  }
  \subsection{Built-in association measures}{
    The names of the following measures can be abbreviated to a unique prefix. Equations are given in the notation of Evert (2008).
    
    \describe{
      \item{\code{frequency} (default)}{
        Co-occurrence \bold{frequency}: \deqn{
          O_{11}
        }{ O11 }
        Use this association measure to operate on raw, unweighted co-occurrence frequency data.
      }
      \item{\code{MI}}{
        \bold{(Pointwise) Mutual Information}, a log-transformed version of the ratio between observed and expected co-occurrence frequency: \deqn{
          \log_2 \frac{O_{11}}{E_{11}}
        }{ log2(O11 / E11) }
        Pointwise MI has a very strong bias towards pairs with low expected co-occurrence frequency (because of \eqn{E_{11}}{E11} in the denominator). It should only be applied if low-frequency targets and features have been removed from the DSM.
        
        The sparse version of MI (with negative scores cut off at 0) is sometimes referred to as "positive pointwise Mutual Information" (\bold{PPMI}) in the literature.
      }
      \item{\code{log-likelihood}}{
        The \eqn{G^2} statistic of a likelihood ratio test for independence of rows and columns in a contingency table, which is very popular in computational linguistics under the name \bold{log-likelihood}: \deqn{
          \pm 2 \left( \sum_{ij} O_{ij}\cdot \log \frac{O_{ij}}{E_{ij}} \right)
        }{ ± 2 * ( SUM[ij] Oij * log(Oij / Eij) ) }
        This implementation computes \emph{signed} association scores, which are negative iff \eqn{O_{11} < E_{11}}{O11 < E11}.
        Log-likelihood has a strong bias towards high co-occurrence frequency and often produces a highly skewed distribution of scores. It may therefore be advisable to combine it with an additional \code{log} transformation.
      }
      \item{\code{simple-ll}}{
        Simple \bold{log-likelihood} (Evert 2008, p. 1225): \deqn{
          \pm 2 \left( O_{11}\cdot \log \frac{O_{11}}{E_{11}} - (O_{11} - E_{11}) \right)
        }{ ± 2 * ( O11 * log(O11 / E11) - (O11 - E11) ) }
        This measure provides a good approximation to the full log-likelihood measure (Evert 2008, p. 1235), but can be computed much more efficiently. It is also very similar to the \bold{local-MI} measure used by several popular DSMs.

        Like \code{log-likelihood}, this measure computes \emph{signed} association scores and has a strong bias towards high co-occurrence frequency.
      }
      \item{\code{t-score}}{
        The \bold{t-score} association measure, which is popular for collocation identification in computational lexicography: \deqn{
          \frac{O_{11} - E_{11}}{\sqrt{O_{11}}}
        }{ (O11 - E11) / sqrt(O11) }
        T-score is known to filter out low-frequency data effectively.  If used as a non-sparse measure, a \dQuote{discounted} version with \eqn{\sqrt(O + 1)} in the denominator is computed.
      }
      \item{\code{chi-squared}}{
        The \eqn{X^2} statistic of Pearson's \bold{chi-squared} test for independence of rows and columns in a contingency table, with Yates's correction applied: \deqn{
          \pm \frac{
            N \bigl( | O_{11}O_{22} - O_{12} O_{21} | - N/2 \bigr)^2
          }{
            R_1 R_2 C_1 C_2
          }
        }{ ± N * (|O12 * O22 - O12 * O21| - N/2)^2 / (R1 * R2 * C1 * C2) }
        This implementation computes \emph{signed} association scores, which are negative iff \eqn{O_{11} < E_{11}}{O11 < E11}.

        The formula above gives a more compact form of Yates's correction than the familiar sum over the four cells of the contingency table.
      }
      \item{\code{z-score}}{
        The \bold{z-score} association measure, based on a normal approximation to the binomial distribution of co-occurrence by chance: \deqn{
          \frac{O_{11} - E_{11}}{\sqrt{E_{11}}}
        }{ (O11 - E11) / sqrt(E11) }
        Z-score has a strong bias towards pairs with low expected co-occurrence frequency (because of \eqn{E_{11}}{E11} in the denominator). It should only be applied if low-frequency targets and features have been removed from the DSM.
      }
      \item{\code{Dice}}{
        The \bold{Dice coefficient} of association, which corresponds to the harmonic mean of the conditional probabilities \eqn{P(\mathrm{feature}|\mathrm{target})}{P(feature | target)} and \eqn{P(\mathrm{target}|\mathrm{feature})}{P(target | feature)}: \deqn{
          \frac{2 O_{11}}{R_1 + C_1}
        }{ 2 O11 / (R1 + C1) }
        Note that Dice is inherently sparse: it preserves zeroes and does not produce negative scores.
      }
    }
    
    The following additional scoring functions can be selected:
    \describe{
      \item{\code{tf.idf}}{
        The \bold{tf-idf} weighting scheme popular in Information Retrieval: \deqn{
          O_{11}\cdot \log \frac{1}{\mathit{df}}
        }{ O11 * log(1 / df) }
        where \eqn{\mathit{df}}{df} is the relative document frequency of the corresponding feature term and should be provided as a variable \code{df} in the model's column information.  Otherwise, it is approximated by the feature's nonzero count \eqn{n_p}{np} (variable \code{nnzero}) divided by the number \eqn{K} of rows in the co-occurrence matrix: \deqn{
          \mathit{df} = \frac{n_p + 1}{K + 1}
        }{ df = (np + 1) / (K + 1) }
        The discounting avoids division-by-zero errors when \eqn{n_p = 0}{np = 0}.
      } 
      \item{\code{reweight}}{
        Apply scale transformation, column scaling and/or row normalization to previously computed feature scores (from \code{model$S}).  This is the only \code{score} that can be used with a DSM that does not contain raw co-occurrence frequency data.
      }
    }
  }
  \subsection{Sparse association scores}{
    If \code{sparse=TRUE}, negative association scores are cut off at 0 in order to (i) ensure that the scored matrix is non-negative and (ii) preserve sparseness.  The implementation assumes that association scores are always \eqn{\leq 0}{≤ 0} for \eqn{O_{11} = 0}{O11 = 0} in this case and only computes scores for nonzero entries in a sparse matrix.  All built-in association measures satisfy this criterion.
    
    Other researchers sometimes refer to such sparse scores as "positive" measures, most notably positive point-wise Mutual Information (PPMI). Since \code{sparse=TRUE} is the default setting, \code{score="MI"} actually computes the PPMI measure.

    Non-sparse association scores can only be computed if \code{negative.ok=TRUE} and will force a dense matrix representation. For this reason, the default is \code{FALSE} for a sparse co-occurrence matrix and \code{TRUE} for a dense one.  A special setting \code{negative.ok="nonzero"} is provided for those who wish to abuse \code{dsm.score} for collocation analysis.  In combination with \code{sparse=FALSE}, it will allow negative score values, but compute them only for the nonzero entries of a sparse co-occurrence matrix.  For a dense co-occurrence matrix, this setting is fully equivalent to \code{negative.ok=TRUE}.
  }
  \subsection{Scale transformations}{
    Association scores can be re-scaled with an isotonic transformation function that preserves sign and ranking of the scores. This is often done in order to de-skew the distribution of scores or as an approximate binarization (presence vs. absence of features).  The following built-in transformations are available:
    \describe{
      \item{\code{none} (default)}{
        A \bold{linear} transformation leaves association scores unchanged. \deqn{
          f(x) = x
        }
      }
      \item{\code{log}}{
        The \bold{logarithmic} transformation has a strong de-skewing effect.  In order to preserve sparseness and sign of association scores, a signed and discounted version has been implemented. \deqn{
          f(x) = \mathop{\mathrm{sgn}}(x) \cdot \log (|x| + 1)
        }{ f(x) = sgn(x) * log(|x| + 1) }
      }
      \item{\code{root}}{
        The \bold{signed square root} transformation has a mild de-skewing effect. \deqn{
          f(x) = \mathop{\mathrm{sgn}}(x) \cdot \sqrt{|x|}
        }{ f(x) = sgn(x) * sqrt(|x|) }
      }
      \item{\code{sigmoid}}{
        The \bold{sigmoid} transformation produces a smooth binarization where negative values saturate at \eqn{-1}, positive values saturate at \eqn{+1} and zeroes remain unchanged. \deqn{
          f(x) = \tanh x
        }{ f(x) = tanh(x) }
      }
    }
  }
  \subsection{User-defined association measures}{
     Instead of the name of a built-in AM, a function implementing a user-defined measure can be passed in the \code{score} argument. This function will be applied to the co-occurrence matrix in batches of approximately \code{batchsize} elements in order to limit the memory overhead incurred. A user-defined AM can be combined with any of the transformations above, and \code{sparse=TRUE} will cut off all negative scores.
     
     The user function can use any of following arguments to access the contingency tables of observed and expected frequencies, following the notation of Evert (2008):
     \describe{
       \item{\code{O}, \code{E}}{observed and expected co-occurrence frequency}
       \item{\code{R1}, \code{R2}, \code{C1}, \code{C2}}{the row and column marginals of the contingency table}
       \item{\code{N}}{sample size}
       \item{\code{f}, \code{f1}, \code{f2}}{the frequency signature of a target-feature pair, a different notation for \eqn{f = O}{f = O}, \eqn{f_1 = R_1}{f1 = R1} and \eqn{f_2 = C_1}{f2 = C1}}
       \item{\code{O11}, \code{O12}, \code{O21}, \code{O22}}{the contingency table of observed frequencies}
       \item{\code{E11}, \code{E12}, \code{E21}, \code{E22}}{the contingency table of expected frequencies}
       \item{\code{rows}}{a data frame containing information about the target items (from the \code{rows} element of \code{model})}
       \item{\code{cols}}{a data frame containing information about the feature items (from the \code{cols} element of \code{model})}
       \item{\code{\dots}}{must be specified to ignore unused arguments}
    }
    Except for \code{rows} and \code{cols}, all these arguments will be numeric vectors of the same lengths or scalar values (\code{N}), and the function must return a numeric vector of the same length.
     
    For example, the built-in Mutual Information measure could also be implemented with the user function
    \preformatted{    my.MI <- function (O, E, ...) log2(O / E) }
    and tf.idf scoring could be implemented as follows, provided that the feature information table \code{model$cols} contains a column \code{df} with relative document frequencies:
    \preformatted{    my.tfidf <- function (O11, cols, ...) O11 * log(1 / cols$df)
    dsm.score(model, score=my.tfidf)}
     
     
     \bold{Warning:} User-defined AMs are much less efficient than the built-in measures and should only be used on large data sets if there is a good reason to do so. Increasing \code{batchsize} may speed up the computation to some degree at the expense of bigger memory overhead.
  }
}
\value{
  Either an updated DSM model of class \code{dsm} (default) or the matrix of (scaled and normalised) association scores (\code{matrix.only=TRUE}).
  
  Note that updating DSM models may require a substantial amount of temporary memory (because of the way memory management is implemented in \R).  This can be problematic when running a 32-bit build of \R or when dealing with very large DSM models, so it may be better to return only the scored matrix in such cases.
}
\references{

  More information about assocation measures and the notation for contingency tables can be found at \url{http://www.collocations.de/} and in

  Evert, Stefan (2008). Corpora and collocations. In A. Lüdeling and M. Kytö (eds.), \emph{Corpus Linguistics. An International Handbook}, chapter 58, pages 1212--1248. Mouton de Gruyter, Berlin, New York.

}
\seealso{

  \code{\link{dsm}}
  
}
\examples{
model <- DSM_TermTerm
model$M # raw co-occurrence matrix
  
model <- dsm.score(model, score="MI")
round(model$S, 3) # PPMI scores
  
model <- dsm.score(model, score="reweight", transform="sigmoid")
round(model$S, 3) # additional sigmoid transformation

## user-defined scoring functions can implement additional measures,
## e.g. the conditional probability Pr(feature | target) as a percentage
my.CP <- function (O11, R1, ...) 100 * O11 / R1  # "..." is mandatory
model <- dsm.score(model, score=my.CP)
round(model$S, 3)

## shifted PPMI (with k = 2) creates all-zero rows and columns
model <- dsm.score(model, score=function (O, E, ...) log2(O / E) - 2,
                   normalize=TRUE, update.nnzero=TRUE)
round(model$S, 3) # normalization preserves all-zero rows
## use subset to remove such rows and columns
m2 <- subset(model, nnzero > 0, nnzero > 0) # must have updated nnzero counts
round(m2$S, 3)

\dontrun{
# visualization of the scale transformations implemented by dsm.score
x <- seq(-2, 4, .025)
plot(x, x, type="l", lwd=2, xaxs="i", yaxs="i", xlab="x", ylab="f(x)")
abline(h=0, lwd=0.5); abline(v=0, lwd=0.5)
lines(x, sign(x) * log(abs(x) + 1), lwd=2, col=2)
lines(x, sign(x) * sqrt(abs(x)), lwd=2, col=3)
lines(x, tanh(x), lwd=2, col=4)
legend("topleft", inset=.05, bg="white", lwd=3, col=1:4,
       legend=c("none", "log", "root", "sigmoid"))
}
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
