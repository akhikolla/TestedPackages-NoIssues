\name{plot.ictest}
\alias{plot.ictest}

\title{
Scatterplot Matrix for a ictest Object
}
\description{
For an object of class ictest, plots either the pairwise scatter plot matrix, or the time series plots of the underlying components. The user can choose if only the components considered interesting or all of them should be plotted.
}
\usage{
\method{plot}{ictest}(x, which = "all", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of class ictest}
  \item{which}{if \code{"all"}, then all components of S in the ictest object are plotted. If \code{"k"}, then only the first k components are plotted,
   where the value of \code{k} is taken from the ictest object. This is only meaningful if \code{k} was at least 2.}
  \item{\dots}{other arguments passed on to \code{\link{pairs}} if the components are a numeric matrix or to \code{\link{plot.ts}}, \code{\link[zoo]{plot.zoo}} or \code{\link[xts]{plot.xts}} if the components are from the corresponding class.}
}

\details{
If the component matrix has the class \code{mts}, \code{xts} or \code{zoo}, then a time series plot will be plotted. Otherwise, the pairwise scatter plot matrix will be plotted.

}

\author{
Klaus Nordhausen
}

\seealso{
\code{\link{ggplot.ictest}, \link{pairs}, \link{plot.ts}, \link[zoo]{plot.zoo}, \link[xts]{plot.xts}}
}
\examples{
n <- 200
X <- cbind(rnorm(n, sd = 2), rnorm(n, sd = 1.5), rnorm(n), rnorm(n), rnorm(n))

TestCov <- PCAasymp(X, k = 2)
plot(TestCov)
plot(TestCov, which = "k")
}


\keyword{hplot }
