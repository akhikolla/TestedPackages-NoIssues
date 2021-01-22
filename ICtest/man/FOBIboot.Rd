\name{FOBIboot}
\alias{FOBIboot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Boostrap-based Testing for the Number of Gaussian Components in ICA Using FOBI 
}
\description{
In independent components analysis (ICA) gaussian components are considered as uninteresting.
The function uses boostrappping tests, based on FOBI, to decide if there are \code{p-k} gaussian components where \code{p} is the dimension of the data.
The function offers two different boostrapping strategies.
}
\usage{
FOBIboot(X, k, n.boot = 200, s.boot = "B1")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{X}{a numeric data matrix with p>1 columns.}
  \item{k}{the number of non-gaussian components under the null.}
  \item{n.boot}{number of bootstrapping samples.}
  \item{s.boot}{bootstrapping strategy to be used. Possible values are \code{"B1"}, \code{"B2"}. See details for further information.}
}
\details{
As in \code{\link{FOBIasymp}} the function jointly diagonalizes the regular covariance and the matrix of fourth moments. Note that in this case the matrix of fourth moments
is not made consistent under the normal model by dividing it by \eqn{p+2}, as for example done by the function \code{\link[ICS]{cov4}} where \eqn{p} denotes the dimension
of the data. Therefore the eigenvalues of this generalized eigenvector-eigenvalue problem which correspond to normally distributed components should be \code{p+2}.
Given eigenvalues \eqn{d_1,...,d_p} the function thus orders the components in descending order according to the values of \eqn{(d_i-(p+2))^2}.

Under the null it is then assumed that the first \code{k} interesting components are mutually independent and non-normal and the last \code{p-k} components are gaussian.

Let \eqn{d_1,...,d_p} be the ordered eigenvalues, \eqn{W} the correspondingly ordered unmixing matrix, \eqn{s_i = W (x_i-MU)}{s_i = W (X_i-MU)} the corresponding
source vectors which give the source matrix \eqn{S} which can be decomposed into \eqn{S_1} and \eqn{S_2} where \eqn{S_1} is the matrix with the \eqn{k} non-gaussian components
and \eqn{S_2} the matrix with the gaussian components (under the null).

The test statistic is then \eqn{T = n \sum_{i=k+1}^p (d_i-(p+2))^2}{T = n sum_(i=k+1)^p (d_i-(p+2))^2}

Two possible bootstrap tests are provided for testing that the last \code{p-k} components are gaussian and independent from the first k components:

\enumerate{
\item \code{s.boot="B1"}: 
The first strategy has the followong steps:
\enumerate{
\item Take a bootstrap sample \eqn{S_1^*} of size \eqn{n} from \eqn{S_1}.
\item Take a bootstrap sample \eqn{S_2^*} consisting of a matrix of standard normally distributed elements.
\item Combine \eqn{S^*=(S_1^*, S_2^*)} and create \eqn{X^*= S^* W}.
\item Compute the test statistic based on \eqn{X^*}. 
\item Repeat the previous steps \code{n.boot} times.
}

Note that in this bootstrapping test the assumption of ''independent components'' is not used, it is only used that the last \eqn{p-k} components are gaussian and independent from the first \eqn{k} components. Therefore this strategy can be applied in an independent component analysis (ICA) framework
and in a non-gaussian components analysis (NGCA) framework.

\item \code{s.boot="B2"}: 
The second strategy has the following steps:
\enumerate{
\item Take a bootstrap sample \eqn{S_1^*} of size \eqn{n} from \eqn{S_1} where the subsampling is done separately for each independent component.
\item Take a bootstrap sample \eqn{S_2^*} consisting of a matrix of standard normally distributed elemenets.
\item Combine \eqn{S^*=(S_1^*, S_2^*)} and create \eqn{X^*= S^* W}.
\item Compute the test statistic based on \eqn{X^*}. 
\item Repeat the previous steps \code{n.boot} times.
}

This bootstrapping strategy assumes a full ICA model and cannot be used in an NGCA framework.
}


}
\value{
A list of class ictest inheriting from class htest containing:
  \item{statistic}{the value of the test statistic.}
  \item{p.value }{the p-value of the test.}
  \item{parameter }{the number of boostrapping samples used to obtain the p-value.}
  \item{method}{character string which test was performed.}
  \item{data.name}{character string giving the name of the data.}
  \item{alternative}{character string specifying the alternative hypothesis.}
  \item{k}{the number or non-gaussian components used in the testing problem.}
  \item{W}{the transformation matrix to the independent components. Also known as unmixing matrix.}
  \item{S}{data matrix with the centered independent components.}
  \item{D}{the underlying FOBI eigenvalues.}
  \item{MU}{the location of the data which was substracted before calculating the independent components.}
  \item{s.boot}{character string which boostrapping strategy was used.}
}

\references{
\cite{Nordhausen, K., Oja, H. and Tyler, D.E. (2016), Asymptotic and bootstrap tests for subspace dimension, <https://arxiv.org/abs/1611.04908>.}

\cite{Nordhausen, K., Oja, H., Tyler, D.E. and Virta, J. (2017), Asymptotic and Bootstrap Tests for the Dimension of the Non-Gaussian Subspace, Signal Processing Letters, 24, 887--891. <doi:10.1109/LSP.2017.2696880 >.}
}
\author{
Klaus Nordhausen
}


\seealso{
\code{\link[JADE]{FOBI}},  \code{\link{FOBIasymp}}
}
\examples{
n <- 1500
S <- cbind(runif(n), rchisq(n, 2), rexp(n), rnorm(n), rnorm(n), rnorm(n))
A <- matrix(rnorm(36), ncol = 6)
X <- S \%*\% t(A)

FOBIboot(X, k = 2)
FOBIboot(X, k = 3, s.boot = "B1")
FOBIboot(X, k = 0, s.boot = "B2")
}

\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
\keyword{ multivariate }% __ONLY ONE__ keyword per line
