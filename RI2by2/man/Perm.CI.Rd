\name{Perm.CI}
\alias{Perm.CI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Permutation test confidence interval for a treatment effect on a
  binary outcome}
\description{
Computes permutation-based confidence intervals for the
average treatment effect on a binary outcome in an experiment where
\eqn{m} of \eqn{n} individuals are randomized to treatment by design.}
\usage{
Perm.CI(data, level, nperm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{observed 2 by 2 table in matrix form where row 1 is
    the treatment assignment Z=1 and column 1 is the binary outcome Y=1}
  \item{level}{significance level of hypothesis tests, i.e., method yields a 100(1-\code{level})\% confidence interval}
  \item{nperm}{number of randomizations to perform for each
    hypothesis test}
}
\details{The permutation confidence interval results from inverting
  \eqn{O(n^{4})} hypothesis tests where \eqn{n} is the total number of
  observations in the observed 2 by 2 table.  For each hypothesis test,
  if \eqn{n \choose m} is less than or equal to \code{nperm}, \eqn{n \choose m}
  randomizations are performed, but if \eqn{n
    \choose m} is greater than \code{nperm}, a random sample with replacement of \code{nperm} randomizations
  are performed.
  }
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
 Perm.CI(ex,0.05,100)
}

\keyword{permutation test}
\keyword{randomization inference}
