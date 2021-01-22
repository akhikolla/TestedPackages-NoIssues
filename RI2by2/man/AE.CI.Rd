\name{AE.CI}
\alias{AE.CI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Attributable effects based confidence interval for a treatment effect on a
  binary outcome}
\description{Computes the attributable effects based confidence interval for the
  average treatment effect on a binary outcome in an experiment where
  \eqn{m} of \eqn{n} individuals are randomized to treatment by design.}
\usage{
AE.CI(data, level)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{observed 2 by 2 table in matrix form where row 1 is
    the treatment assignment Z=1 and column 1 is the binary outcome Y=1}
  \item{level}{significance level of hypothesis tests, i.e., method yields a 100(1-\code{level})\% confidence interval}
}
\details{The attributable effects based confidence interval from inverting \eqn{n+2} hypothesis tests.}
\value{
\item{tau.hat}{estimated average treatment effect}
\item{lower}{lower bound of confidence interval}
\item{upper}{upper bound of confidence interval}
}
\references{
Rigdon, J.R. and Hudgens, M.G. (2015).
Randomization inference for treatment effects on a binary outcome.
\emph{Statistics in Medicine, 34(6), 924-935}.}

\author{Joseph Rigdon \email{jrigdon@stanford.edu}}

\examples{
 ex = matrix(c(8,2,3,7),2,2,byrow=TRUE)
 AE.CI(ex,0.05)
}

\keyword{attributable effects}
\keyword{randomization inference}
