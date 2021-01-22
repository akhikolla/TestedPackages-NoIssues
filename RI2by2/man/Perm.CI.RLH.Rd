\name{Perm.CI.RLH}
\alias{Perm.CI.RLH}
\title{Permutation test confidence interval for a treatment effect on a binary outcome}
\description{
Computes permutation-based confidence intervals for the
average treatment effect on a binary outcome in an experiment where
\eqn{m} of \eqn{n} individuals are randomized to treatment by design.
This function is based on the modified approach (RLH) in
Rigdon, Loh and Hudgens (forthcoming).
The Chiba (2015) and Blaker (2000) intervals are also returned.
There is an additional option of specifying the maximum number of
hypothesis tests to be carried out.
}
\usage{
Perm.CI.RLH(data, level, verbose=FALSE, total_tests=NA)
}
\arguments{
\item{data}{observed 2 by 2 table in matrix form where row 1 is
    the treatment assignment Z=1 and column 1 is the binary outcome Y=1}
\item{level}{significance level of hypothesis tests, i.e., method yields a 100(1-\code{level})\% confidence interval}

\item{verbose}{If \code{TRUE}, returns an additional data frame listing
all the values of
\eqn{(n_{11},n_{10},n_{01},n_{00})}{(n11,n10,n01,n00)}
tested, and the corresponding p-values; default = \code{FALSE}.}

\item{total_tests}{maximum number of hypotheses to be tested in total,
with a minimum of two for each possible value of
\eqn{(n_{10}-n_{01})/n;}{(n10-n01)/n;}
default = \code{NA}.
By default, all hypotheses are evaluated
until the minimum and maximum values of
\eqn{(n_{10}-n_{01})/n}{(n10-n01)/n}
with p-values \eqn{\ge} \code{level} (or \code{level}/2 for the Chiba intervals)
are found.}

}
\value{A list with the following items:
\item{Chiba}{Chiba confidence interval}
\item{RLH}{RLH confidence interval}
\item{Blaker}{Blaker confidence interval}
\item{tau.hat}{estimated average treatment effect}
\item{p_values}{if \code{verbose=TRUE}, a data frame with all the p-values from the hypothesis tests; default=\code{FALSE}}
}
\references{
Rigdon, J.R. and Hudgens, M.G. (2015).
Randomization inference for treatment effects on a binary outcome.
\emph{Statistics in Medicine, 34(6), 924-935}.

Chiba, Y. (2015).
Exact tests for the weak causal null hypothesis on a binary outcome in randomized trials.
\emph{Journal of Biometrics & Biostatistics, 6(244)}.

Chiba, Y. (2016).
A note on exact confidence interval for causal effects on a binary outcome in randomized trials.
\emph{Statistics in Medicine, 35(10), 1739-1741}.

Blaker, H. (2000).
Confidence curves and improved exact confidence intervals for discrete distributions.
\emph{Canadian Journal of Statistics, 28(4), 783-798}.

Rigdon, J.R., Loh W.W. and Hudgens, M.G. (forthcoming).
Response to comment on "Randomization inference for treatment effects on a binary outcome."
}

\author{Wen Wei Loh \email{wwloh@email.unc.edu}}

\examples{
ex = matrix(c(11,1,7,21),2,2,byrow=TRUE)
Perm.CI.RLH(ex,0.05)

ex = matrix(c(7,5,1,27),2,2,byrow=TRUE)
Perm.CI.RLH(ex,0.05)
Perm.CI.RLH(ex,0.05, verbose=TRUE)

ex = matrix(c(33,15,11,37),2,2,byrow=TRUE)
Perm.CI.RLH(ex,0.05, total_tests=1000)
\donttest{Perm.CI.RLH(ex,0.05)

ex = matrix(c(40,60,15,85),2,2,byrow=TRUE)
Perm.CI.RLH(ex,0.05, total_tests=1000)}
\dontrun{Perm.CI.RLH(ex,0.05)}
}