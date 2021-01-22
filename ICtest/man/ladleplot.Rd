\name{ladleplot}
\alias{ladleplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ladle Plot for an Object of Class ladle 
}
\description{
The ladle plot is a measure to decide about the number of interesting components. Of interest for the ladle criterion is the minimum.
The function here offers however also to plot other criterion values which are part of the actual ladle criterion. 
}
\usage{
ladleplot(x, crit = "gn", type="l", ylab = crit, 
          xlab = "component", main = deparse(substitute(x)), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class ladle.}
  \item{crit}{the criterion to be plotted, options are \code{"gn"}, \code{"fn"}, \code{"phin"} and \code{"lambda"}.}
  \item{type}{plotting type.}
  \item{ylab}{default ylab value.}
  \item{xlab}{default xlab value.}
  \item{main}{default title.}
  \item{\dots}{other arguments for the plotting functions.}
}
\details{
The main criterion of the ladle is the scaled sum of the eigenvalues and the measure of variation of the eigenvectors up to the component of interest.

The sum is denoted \code{"gn"} and the individual parts are \code{"fn"} for the measure of the eigenvector variation and \code{"phin"} for the scaled eigenvalues.
The last option \code{"lambda"} corresponds to the unscaled eigenvalues yielding then a screeplot.
}

\references{
\cite{Luo, W. and Li, B. (2016), Combining Eigenvalues and Variation of Eigenvectors for Order Determination, Biometrika, 103. 875--887. <doi:10.1093/biomet/asw051>}
}
\author{
Klaus Nordhausen
}

\seealso{
\code{\link{FOBIladle}}, \code{\link{PCAladle}}, \code{\link{SIRladle}}
}
\examples{
n <- 1000
X <- cbind(rexp(n), rt(n,5), rnorm(n), rnorm(n), rnorm(n), rnorm(n))
test <- FOBIladle(X)
ladleplot(test)
ladleplot(test, crit="fn")
ladleplot(test, crit="phin")
ladleplot(test, crit="lambda")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }% use one of  RShowDoc("KEYWORDS")
