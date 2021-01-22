\name{SIRasymp}
\alias{SIRasymp}

\title{
Testing the Subspace Dimension for Sliced Inverse Regression.
}
\description{
Using the two scatter matrices approach (SICS) for sliced inversion regression (SIR), the function tests
if the last \code{p-k} components have zero eigenvalues, where \code{p} is the number of explaining variables. Hence the assumption is that the first \code{k} 
components are relevant for modelling the response \code{y} and the remaining components are not.
}
\usage{
SIRasymp(X, y, k, h = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrix of explaining variables.}
  \item{y}{a numeric vector specifying the response.}
  \item{k}{the number of relevant components under the null hypothesis.}
  \item{h}{the number of slices used in SIR. Passed on to function \code{\link{covSIR}}.}
  \item{\dots}{other arguments passed on to \code{\link{covSIR}}.}
}
\details{
Under the null the first \code{k} eigenvalues contained in \code{D} are non-zero and the remaining \code{p-k} are zero.

For a sample of size \eqn{n}, the test statistic \eqn{T} is then n times the sum of these last p-k eigenvalue and has under the null a chisquare distribution with \eqn{(p-k)(h-k-1)} degrees of freedom,
therefore it is required that \eqn{k < h-1}.
}
\value{
A list of class ictest inheriting from class htest containing:
  \item{statistic}{the value of the test statistic.}
  \item{p.value }{the p-value of the test.}
  \item{parameter }{the degrees of freedom of the test.}
  \item{method}{character string which test was performed.}
  \item{data.name}{character string giving the name of the data.}
  \item{alternative}{character string specifying the alternative hypothesis.}
  \item{k}{the number of non-zero eigenvalues used in the testing problem.}
  \item{W}{the transformation matrix to the underlying components.}
  \item{S}{data matrix with the centered underlying components.}
  \item{D}{the underlying eigenvalues.}
  \item{MU}{the location of the data which was substracted before calculating the components.}
}
\references{
\cite{Nordhausen, K., Oja, H. and Tyler, D.E. (2016), Asymptotic and bootstrap tests for subspace dimension, <https://arxiv.org/abs/1611.04908>.}
}
\author{
Klaus Nordhausen
}



\seealso{
\code{\link{covSIR}},  \code{\link{SIRboot}}
}
\examples{
X <- matrix(rnorm(1000), ncol = 5)
eps <- rnorm(200, sd = 0.1)
y <- 2 + 0.5 * X[, 1] + 2 * X[, 3] + eps
  
SIRasymp(X, y, k = 0) 
SIRasymp(X, y, k = 1)    
}

\keyword{ multivariate }
\keyword{ htest }
