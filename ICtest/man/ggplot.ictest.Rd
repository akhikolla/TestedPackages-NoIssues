\name{ggplot.ictest}
\alias{ggplot.ictest}

\title{
Scatterplot Matrix for a ictest Object using ggplot2
}
\description{
For an object of class ictest, plots either the pairwise scatter plot matrix using \code{ggpairs} from GGally, or the time series plots of the underlying components using ggplot2. The user can choose if only the components considered interesting or all of them should be plotted. Aesthetics can be passed to ggpairs as well.
}
\usage{
\method{ggplot}{ictest}(data, mapping = aes(), mapvar = NULL, which = "all", ..., 
       environment=parent.frame())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{object of class ictest}
  \item{mapping}{aesthetic mapping, see documentation for \code{\link{ggpairs}}. If \code{x} has the class \code{mts} then this argument is not used.}
  \item{mapvar}{data.frame of the external variables used by the aesthetic mappings. If \code{x} has the class \code{mts} then this argument is not used.}
  \item{which}{if \code{"all"}, then all components of S in the ictest object are plotted. If \code{"k"}, then only the first k components are plotted,
   where the value of \code{k} is taken from the ictest object. This is only meaningful if \code{k} was at least 2.}
  \item{...}{arguments passed on to \code{\link{ggpairs}}. If the component matrix has the class \code{mts}, \code{xts} or \code{zoo} then this argument is not used.}
  \item{environment}{not used but needed for consistency.}
}

\details{
If the component matrix has the class \code{mts}, \code{xts} or \code{zoo} then a time series plot will be plotted using ggplot2. Otherwise, a pairwise scatter plot matrix will be plotted using GGally.

}


\author{
Klaus Nordhausen, Joni Virta
}

\seealso{
\code{\link{plot.ictest}, \link{pairs}}
}
\examples{
# The data
X <- iris[, 1:4]

# The aesthetics variables
mapvar <- data.frame(iris[, 5])
colnames(mapvar) <- "species"

TestCov <- PCAasymp(X, k = 2)
ggplot(TestCov)
ggplot(TestCov, aes(color = species), mapvar = mapvar, which = "k")
}


\keyword{hplot }
