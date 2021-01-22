\name{plot.eval.similarity.correlation}
\alias{plot.eval.similarity.correlation}
\alias{print.eval.similarity.correlation}
\title{
  Printing and Plotting Similarity Correlation Evaluation Results (wordspace)
}
\description{

  Suitable printing and visualization of evaluation results from \code{\link{eval.similarity.correlation}}.
  The \code{print} method displays an evaluation summary (stored in attribute \code{"eval.result"}) after the full data frame.
  The \code{plot} method displays a scatterplot of gold standard ratings against DSM distances with optional
  regression line (\code{\link{lowess}}), a summary of evaluation results at the top, and various other formatting options.

}
\usage{
  \method{print}{eval.similarity.correlation}(x, \dots)
  
  \method{plot}{eval.similarity.correlation}(x, y, line = TRUE,
     categories = NULL, cat.col = NA, cat.legend = "bottomleft",
     pch = 20, cex = 1, xlim = NULL, ylim = NULL,
     xlab = "human rating", ylab = "distributional model",
     main = attr(x, "taskname"), \dots)
}
\arguments{

  \item{x}{detailed evaluation report from \code{\link{eval.similarity.correlation}} (with \code{details=TRUE})}
  
  \item{y}{unused, must not be specified}
  
  \item{line}{if \code{TRUE}, a non-linear regression line is added to the plot in order to indicate the precise relationship between gold standard ratings and DSM distances (based on the \code{\link{lowess}} smoother)}

  \item{categories}{a factor with one entry for each word pair in the evaluation task.
    If specified, points in the scatterplot are colour-coded according to the categories of the corresponding word pairs.
    Note that \code{categories} is evaluated within the data frame \code{x}, so any column of \code{x}
    can directly be used as a variable.
  }

  \item{cat.col}{a vector of colours to be used for the different categories (defaults to the standard palette built into \R)}

  \item{cat.legend}{corner of the plot in which to display the legend box for category colours}

  \item{pch, cex}{symbol used for plotting and its size}
  
  \item{xlim, ylim}{range of values shown on the x-axis and y-axis, respectively}
  
  \item{xlab, ylab, main}{axis labels and main title of the plot}

  \item{\dots}{all other arguments are passed to the scatterplot function (\code{\link{plot.default}}) and can be used to set additional graphics parameters}

}
\details{

  Word pairs not found in the DSM are always shown as empty boxes in the scatterplot, regardless of the \code{pch} parameter.

}
\seealso{
  \code{\link{eval.similarity.correlation}}
}
\examples{
\dontrun{
plot(eval.similarity.correlation(WordSim353, DSM_Vectors, details=TRUE))
}
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{hplot}
