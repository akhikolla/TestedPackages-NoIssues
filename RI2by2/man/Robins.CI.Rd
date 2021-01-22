\name{Robins.CI}
\alias{Robins.CI}
\title{Asymptotic confidence interval for a treatment effect on a
  binary outcome}
\description{Computes the Robins (1988) confidence interval for the
  average treatment effect on a binary outcome in an experiment where
  \eqn{m} of \eqn{n} individuals are randomized to treatment by design.}
\usage{
Robins.CI(data, level)
}

\arguments{
  \item{data}{observed 2 by 2 table in matrix form where row 1 is
    the treatment assignment Z=1 and column 1 is the binary outcome Y=1}
  \item{level}{significance level of hypothesis tests, i.e., method yields a 100(1-\code{level})\% confidence interval}
}

\details{The Robins (1988) confidence interval is similar in form to the well known Wald confidence interval for a difference in proportions, but is guaranteed to have smaller width.}
\value{
\item{tau.hat}{estimated average treatment effect}
\item{lower}{lower bound of confidence interval}
\item{upper}{upper bound of confidence interval}
}

\references{
Robins, J.M. (1988). Confidence intervals for causal parameters.
\emph{Statistics in Medicine, 7(7), 773-785}.}

\author{Joseph Rigdon \email{jrigdon@stanford.edu}}

\examples{
#Example 1 from Robins (1988)
ex = matrix(c(40,60,15,85),2,2,byrow=TRUE)
Robins.CI(ex,0.05)
}

\keyword{randomization inference}

