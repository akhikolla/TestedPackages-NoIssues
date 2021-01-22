\name{PCAasymp}
\alias{PCAasymp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Testing for Subsphericity using the Covariance Matrix or Tyler's Shape Matrix
}
\description{
The function tests, assuming an elliptical model, that the last \code{p-k} eigenvalues of
a scatter matrix are equal and the \code{k} interesting components are those with a larger variance. 
The scatter matrices that can be used here are the regular covariance matrix and Tyler's shape matrix.
}
\usage{
PCAasymp(X, k, scatter = "cov", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrix with p>1 columns.}
  \item{k}{the number of eigenvalues larger than the equal ones. Can be between 0 and p-2.}
  \item{scatter}{the scatter matrix to be used. Can be \code{"cov"} or \code{"tyler"}. For \code{"cov"} the regular covariance matrix is computed and for
  \code{"tyler"} the function \code{\link[ICSNP]{HR.Mest}} is used to compute Tyler's shape matrix.}
  \item{\dots}{arguments passed on to \code{\link[ICSNP]{HR.Mest}} if \code{scatter = "tyler"}.}
}


\details{
The functions assumes an elliptical model and tests if the last \eqn{p-k} eigenvalues of PCA are equal. PCA can here be either be based on the regular covariance matrix or on Tyler's shape matrix.

For a sample of size \eqn{n}, the test statistic is
\deqn{T = n / (2 \bar{d}^2 \sigma_1) \sum_{k+1}^p (d_i - \bar{d})^2,}{T = n / (2 d.bar^2 sigma_1) \sum_(k+1)^p (d_i - d.bar)^2,} 
where \eqn{\bar{d}}{d.bar} is the mean of the last \eqn{p-k} PCA eigenvalues. 

The constant \eqn{\sigma_1}{sigma_1} is for the regular covariance matrix estimated from the data whereas for Tyler's shape matrix it is simply a function of the dimension of the data.

The test statistic has a limiting chisquare distribution with \eqn{(p-k-1)(p-k+2)/2} degrees of freedom.

Note that the regular covariance matrix is here divided by \eqn{n} and not by \eqn{n-1}. 
}

\value{
A list of class ictest inheriting from class htest containing:
  \item{statistic}{the value of the test statistic.}
  \item{p.value }{the p-value of the test.}
  \item{parameter }{the degrees of freedom of the test.}
  \item{method}{character string which test was performed.}
  \item{data.name}{character string giving the name of the data.}
  \item{alternative}{character string specifying the alternative hypothesis.}
  \item{k}{the number or larger eigenvalues used in the testing problem.}
  \item{W}{the transformation matrix to the principal components.}
  \item{S}{data matrix with the centered principal components.}
  \item{D}{the underlying eigenvalues.}
  \item{MU}{the location of the data which was substracted before calculating the principal components.}
  \item{SCATTER}{the computed scatter matrix.}
  \item{sigma1}{the asymptotic constant needed for the asymptotic test.}
}

\references{
\cite{Nordhausen, K., Oja, H. and Tyler, D.E. (2016), Asymptotic and bootstrap tests for subspace dimension, <https://arxiv.org/abs/1611.04908>.}
}
\author{
Klaus Nordhausen
}

\seealso{
\code{\link[ICSNP]{HR.Mest}}, \code{\link{PCAboot}}
}
\examples{
n <- 200
X <- cbind(rnorm(n, sd = 2), rnorm(n, sd = 1.5), rnorm(n), rnorm(n), rnorm(n))

TestCov <- PCAasymp(X, k = 2)
TestCov
TestTyler <- PCAasymp(X, k = 1, scatter = "tyler")
TestTyler
}

\keyword{ htest }
\keyword{ multivariate }
