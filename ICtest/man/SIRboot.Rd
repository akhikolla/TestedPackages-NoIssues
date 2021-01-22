\name{SIRboot}
\alias{SIRboot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Testing the Subspace Dimension for Sliced Inverse Regression Using Bootstrapping.
}
\description{
Using the two scatter matrices approach (SICS) for sliced inversion regression (SIR) the function tests
if the last \code{p-k} components have zero eigenvalues, where \code{p} is the number of explaining variables. Hence the assumption is that the first \code{k} 
components are relevant for modelling the response \code{y} and the remaining components are not. The function performs bootstrapping to obtain a p-value.
}
\usage{
SIRboot(X, y, k, h = 10, n.boot = 200, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrix of explaining variables.}
  \item{y}{a numeric vector specifying the response.}
  \item{k}{the number of relevant components under the null hypothesis.}
  \item{h}{the number of slices used in SIR. Passed on to function \code{\link{covSIR}}.}
  \item{n.boot}{number of bootstrapping samples.}
  \item{\dots}{other arguments passed on to \code{\link{covSIR}}.}
}
\details{
Under the null hypthesis the last p-k eigenvalue as given in D are zero. The test statistic is then the sum of these eigenvalues.

Denote W as the transformation matrix to the supervised invariant coordinates (SIC) \eqn{s_i}, \eqn{i=1,\ldots,n}, i.e.
\deqn{s_i = W (x_i-MU),}{s_i = W (X_i-MU),}
where \code{MU} is the location.

Let \eqn{S_1}{S_1} be the submatrix of the SICs which are relevant and \eqn{S_2}{S_2} the submatrix of the SICs which are irrelevant for the response y under the null.

The boostrapping has then the following steps:

\enumerate{
\item Take a boostrap sample \eqn{(y^*, S_1^*)} of size \eqn{n} from \eqn{(y, S_1)}.
\item Take a boostrap sample \eqn{S_2^*} of size \eqn{n} from \eqn{S_2}.
\item Combine \eqn{S^*=(S_1^*, S_2^*)} and create \eqn{X^*= S^* W}.
\item Compute the test statistic based on \eqn{X^*}. 
\item Repeat the previous steps \code{n.boot} times.
}

}

\value{
A list of class ictest inheriting from class htest containing:
  \item{statistic}{the value of the test statistic.}
  \item{p.value }{the p-value of the test.}
  \item{parameter }{the number of boostrapping samples used to compute the p-value.}
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
\code{\link{covSIR}},  \code{\link{SIRasymp}}
}
\examples{
X <- matrix(rnorm(1000), ncol = 5)
eps <- rnorm(200, sd = 0.1)
y <- 2 + 0.5 * X[, 1] + 2 * X[, 3] + eps
  
SIRboot(X, y, k = 0) 
SIRboot(X, y, k = 1)    
}

\keyword{ multivariate }
\keyword{ htest }
