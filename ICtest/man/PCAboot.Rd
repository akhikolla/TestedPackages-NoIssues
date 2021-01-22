\name{PCAboot}
\alias{PCAboot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap-Based Testing for Subsphericity 
}
\description{
The function tests, assuming an elliptical model, that the last \code{p-k} eigenvalues of
a scatter matrix are equal and the \code{k} interesting components are those with a larger variance. 
To obtain p-values two different bootstrapping strategies are available and the user can provide the scatter matrix to be used
as a function.
}
\usage{
PCAboot(X, k, n.boot = 200, s.boot = "B1", S = MeanCov, Sargs = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrix with p>1 columns.}
  \item{k}{the number of eigenvalues larger than the equal ones. Can be between 0 and p-2.}
  \item{n.boot}{number of bootstrapping samples.}
  \item{s.boot}{bootstrapping strategy to be used. Possible values are \code{"B1"}, \code{"B2"}. See details for further information.}
  \item{S}{A function which returns a list that has as its first element a location vector and as the second element the scatter matrix.}
  \item{Sargs}{list of further arguments passed on to the function specified in \code{S}.}
}

\details{
Here the function \code{S} needs to return a list where the first argument is a location vector and the second one a scatter matrix.

The location is used to center the data and the scatter matrix is used to perform PCA.

Consider X as the centered data and denote by W the transformation matrix to the principal components. The corresponding eigenvalues
from PCA are \eqn{d_1,...,d_p}. Under the null, \eqn{d_k > d_{k+1} = ... = d_{p}}{d_k > d_(k+1) = ... = d_p}.
Denote further by \eqn{\bar{d}}{d.bar} the mean of the last \code{p-k} eigenvalues and by \eqn{D^* = diag(d_1,...,d_k,\bar{d},...,\bar{d})}{D* = diag(d_1,...,d_k,d.bar,...,d.bar)} a \eqn{p \times p}{p x p} diagonal matrix. Assume that \eqn{S} is the matrix with principal components which can be decomposed into \eqn{S_1} and \eqn{S_2} where 
\eqn{S_1} contains the k interesting principal components and \eqn{S_2} the last \eqn{p-k} principal components.

For a sample of size \eqn{n}, the test statistic used for the boostrapping tests is
\deqn{T = n / (\bar{d}^2) \sum_{k+1}^p (d_i - \bar{d})^2.}{T = n / (d.bar^2) sum_(k+1)^p (d_i - d.bar)^2.} 

The function offers then two boostrapping strategies:

\enumerate{
\item \code{s.boot="B1"}: 
The first strategy has the following steps:
\enumerate{
\item Take a bootstrap sample \eqn{S^*} of size \eqn{n} from \eqn{S} and decompose it into \eqn{S_1^*} and \eqn{S_2^*}.
\item Every observation in \eqn{S_2^*} is transformed with a different random orthogonal matrix.
\item Recombine \eqn{S^*=(S_1^*, S_2^*)} and create \eqn{X^*= S^* W}.
\item Compute the test statistic based on \eqn{X^*}. 
\item Repeat the previous steps \code{n.boot} times.
}
\item \code{s.boot="B2"}: 
The second strategy has the following steps:
\enumerate{
\item Scale each principal component using the matrix \eqn{D}, i.e. \eqn{Z = S D}. 
\item Take a bootstrap sample \eqn{Z^*} of size \eqn{n} from \eqn{Z}.
\item Every observation in \eqn{Z^*} is transformed with a different random orthogonal matrix.
\item Recreate \eqn{X^*= Z^* {D^*}^{-1} W}.
\item Compute the test statistic based on \eqn{X^*}. 
\item Repeat the previous steps \code{n.boot} times.
}
To create the random orthogonal matrices the function \code{\link{rorth}} is used.
}
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
  \item{SCATTER}{The computed scatter matrix.}
  \item{scatter}{character string denoting which scatter function was used.}
  \item{s.boot}{character string denoting which bootstrapping test version was used.}
}

\references{
\cite{Nordhausen, K., Oja, H. and Tyler, D.E. (2016), Asymptotic and bootstrap tests for subspace dimension, <https://arxiv.org/abs/1611.04908>.}
}
\author{
Klaus Nordhausen
}

\seealso{
\code{\link{cov}}, \code{\link[ICS]{MeanCov}}, \code{\link{PCAasymp}}
}
\examples{
n <- 200
X <- cbind(rnorm(n, sd = 2), rnorm(n, sd = 1.5), rnorm(n), rnorm(n), rnorm(n))

TestCov <- PCAboot(X, k = 2)
TestCov

# for demonstration purpose the n.boot is chosen small, should be larger in real applications
TestTM <- PCAboot(X, k = 1, n.boot=80, s.boot = "B2", S = "tM", Sargs = list(df=2))
TestTM

}

\keyword{ htest }
\keyword{ multivariate }
