\name{plot.ladle}
\alias{plot.ladle}

\title{
Plotting an Object of Class ladle
}
\description{
An object of class ladle contains always the source components as estimated by the corresponding statistical method.
This function either plots all of the components or only this considered interesting according to the ladle estimate.
}
\usage{
\method{plot}{ladle}(x, which = "all", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class ladle.}
  \item{which}{if \code{"all"}, then all components of S in the ladle object are
          plotted. If \code{"k"}, then only the  k components are
          plotted, which are considered interesting according to the ladle estimator. 
          This is only meaningful if the estimated 'k' is at least 2.}
  \item{\dots}{other arguments passed on to \code{\link{pairs}} if the components are a numeric matrix or to \code{\link{plot.ts}}, \code{\link[zoo]{plot.zoo}} or \code{\link[xts]{plot.xts}} if the components are from the corresponding class.}
}

\details{
If the component matrix has the class \code{mts}, \code{xts} or \code{zoo}, then a time series plot will be plotted. Otherwise, the pairwise scatter plot matrix will be plotted.

}

\author{
Klaus Nordhausen
}

\seealso{
\code{\link{pairs}}, \code{\link{plot.ts}}, \code{\link[zoo]{plot.zoo}}, \code{\link[xts]{plot.xts}}
}
\examples{
n <- 1000
X <- cbind(rexp(n), rt(n,5), rnorm(n), rnorm(n), rnorm(n), rnorm(n))
test <- FOBIladle(X)
plot(test)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }% use one of  RShowDoc("KEYWORDS")
