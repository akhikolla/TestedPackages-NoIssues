\name{PCAschott}
\alias{PCAschott}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Testing for Subsphericity using the Schott's test
}
\description{
The test tests the equality of the last eigenvalues assuming normal distributed data using the regular covariance matrix.
}
\usage{
PCAschott(X, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrix with p>1 columns.}
  \item{k}{the number of eigenvalues larger than the equal ones. Can be between 0 and p-2.}
}
\details{
The functions assumes multivariate normal data and tests if the last \eqn{p-k} eigenvalues of PCA are equal.

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
  \item{MU}{the mean vector of the data which was substracted before calculating the principal components.}
  \item{SCATTER}{the computed covariance matrix matrix.}
}
\references{
\cite{Schott, J.R. (2006), A High-Dimensional Test for the Equality of the Smallest Eigenvalues of a Covariance Matrix, Journal of Multivariate Analysis, 97, 827--843. <doi:10.1016/j.jmva.2005.05.003>}
}
\author{
Klaus Nordhausen
}

\seealso{
\code{\link{PCAasymp}}, \code{\link{PCAboot}}
}
\examples{
n <- 200
X <- cbind(rnorm(n, sd = 2), rnorm(n, sd = 1.5), rnorm(n), rnorm(n), rnorm(n))
PCAschott(X, 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ multivariate }
