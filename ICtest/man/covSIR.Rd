\name{covSIR}
\alias{covSIR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Supervised Scatter Matrix as Used in Sliced Inverse Regression
}
\description{
Sliced Inverse Regression (SIR) can be seen as  special case of Supervised ICS (SICS)
and this function gives the supervised scatter matrix for SIR
}
\usage{
covSIR(X, y, h = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric data matrix.
}
  \item{y}{
a numeric response vector.
}
  \item{h}{
the number of slices.
}
  \item{\dots}{
arguments passed on to \code{\link{quantile}}.
}
}
\details{
This supervised scatter matrix is usually used as the second scatter matrix in SICS to obtain a SIR type supervised linear dimension reduction.
For that purpose \code{covSIR} first divides the response \code{y} into \code{h} slices using the corresponding quantiles as cut points.
Then for each slice the mean vector of \code{X} is computed and the resulting supervised scatter matrix consist of the covariance matrix of these mean vectors.

The function might have problems if the sample size is too small.
}
\value{
a supervised scatter matrix
}
\references{
\cite{Liski, E., Nordhausen, K. and Oja, H. (2014), Supervised invariant coordinate selection, \emph{Statistics: A Journal of Theoretical and Applied Statistics}, \bold{48},  711--731. <doi:10.1080/02331888.2013.800067>.}
 
\cite{Nordhausen, K., Oja, H. and Tyler, D.E. (2016), Asymptotic and bootstrap tests for subspace dimension, <https://arxiv.org/abs/1611.04908>.}
}
\author{
Klaus Nordhausen
}

\seealso{
\code{\link{ics}}
}
\examples{
X <- matrix(rnorm(1000), ncol = 5)
eps <- rnorm(200, sd = 0.1)
y <- 2 + 0.5 * X[, 1] + 2 * X[, 3] + eps

covSIR(X, y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate}% use one of  RShowDoc("KEYWORDS")
