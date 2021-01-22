\name{summary.ladle}
\alias{summary.ladle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing an Object of Class ladle
}
\description{
Summarizes an ladle object
}
\usage{
\method{summary}{ladle}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class ladle.}
  \item{\dots}{further arguments to be passed to or from methods.}
}


\author{
Klaus Nordhausen
}

\seealso{
\code{\link{print.ladle}}, \code{\link{plot.ladle}}, \code{\link{ladleplot}}, \code{\link{FOBIladle}}, \code{\link{PCAladle}}, \code{\link{SIRladle}}
}

\examples{
n <- 1000
X <- cbind(rexp(n), rt(n,5), rnorm(n), rnorm(n), rnorm(n), rnorm(n))

test <- FOBIladle(X)
summary(test)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }% use one of  RShowDoc("KEYWORDS")
\keyword{print }% __ONLY ONE__ keyword per line
