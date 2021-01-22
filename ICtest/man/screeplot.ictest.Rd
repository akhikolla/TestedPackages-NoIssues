\name{screeplot.ictest}
\alias{screeplot.ictest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Screeplot for an ictest Object
}
\description{
Plots the criterion values of an \code{ictest} object against its index number. Two versions of this screeplot are available.
}
\usage{
\method{screeplot}{ictest}(x, type = "barplot", main = deparse(substitute(x)), 
  ylab = "criterion", xlab = "component", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ object of class \code{ictest}.}
  \item{type}{ \code{barplot} if a barplot or \code{lines} if a line plot is preferred.}
  \item{main}{ main title of the plot. }
  \item{ylab}{ y-axis label. }
  \item{xlab}{ x-axis label. }
  \item{\dots}{other arguments for the plotting functions. }  
}


\author{
Klaus Nordhausen
}

\seealso{
\code{\link{ggscreeplot}}
}

\examples{
n <- 200
X <- cbind(rnorm(n, sd = 2), rnorm(n, sd = 1.5), rnorm(n), rnorm(n), rnorm(n))

TestCov <- PCAasymp(X, k = 2)
screeplot(TestCov)
}
\keyword{ hplot }
