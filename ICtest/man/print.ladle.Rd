\name{print.ladle}
\alias{print.ladle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Printing an Object of Class ladle
}
\description{
Basic printing of an object of class ladle. Prints basically everything but the estimated components.
}
\usage{
\method{print}{ladle}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class ladle.}
  \item{\dots}{further arguments to be passed to or from methods.}
}


\author{
Klaus Nordhausen
}
\seealso{
\code{\link{summary.ladle}}, \code{\link{plot.ladle}}, \code{\link{ladleplot}}, \code{\link{FOBIladle}}, \code{\link{PCAladle}}, \code{\link{SIRladle}}
}


\examples{
n <- 1000
X <- cbind(rexp(n), rt(n,5), rnorm(n), rnorm(n), rnorm(n), rnorm(n))
test <- FOBIladle(X)
test
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }% use one of  RShowDoc("KEYWORDS")
\keyword{ print }% __ONLY ONE__ keyword per line
