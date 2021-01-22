\name{NGPP}
\alias{NGPP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Non-Gaussian Projection Pursuit 
}
\description{
Estimates \eqn{k} non-Gaussian signal components using projection pursuit. The projection index can be chosen among convex combinations of squares of one or several standard projection indices used in ICA.
}
\usage{
NGPP(X, k, nl = c("skew", "pow3"), alpha = 0.8, method = "symm", eps = 1e-6,
     verbose = FALSE, maxiter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Numeric matrix with n rows corresponding to the observations and p columns corresponding to the variables.}
  \item{k}{Number of components to estimate, \eqn{1 \leq k \leq p}{1 <= k <= p}.}
  \item{nl}{Vector of non-linearities, a convex combination of the corresponding squared objective functions of which is then used as the projection index. The choices include \code{"skew"} (skewness), \code{"pow3"} (excess kurtosis), \code{"tanh"} (\eqn{log(cosh)}) and \code{"gauss"} (Gaussian function).}
  \item{alpha}{Vector of positive weights between 0 and 1 given to the non-linearities. The length of \code{alpha} should be either one less than the number of non-linearities in which case the missing weight is chosen so that \code{alpha} sums to one, or equal to the number of non-linearities in which case the weights are used as such. No boundary checks for the weights are done.}
  \item{method}{If \code{"symm"} the \code{k} signals are estimated simultaneously (symmetric projection pursuit) and if \code{"defl"} they are estimated one-by-one (deflation-based projection pursuit).}
  \item{eps}{Convergence tolerance.}
  \item{verbose}{If \code{TRUE} the numbers of iterations will be printed.}
  \item{maxiter}{Maximum number of iterations.}
}
\details{
It is assumed that the data is a random sample from the model \eqn{x = m + A s} where the latent vector \eqn{s = (s_1^T, s_2^T)^T}{s = (s_1',s_2')'} consists of \eqn{k}-dimensional non-Gaussian subvector (the signal) and \eqn{p - k}-dimensional Gaussian subvector (the noise) and the components of \eqn{s} are mutually independent. Without loss of generality we further assume that the components of \eqn{s} have zero means and unit variances.

The objective is to estimate an inverse for the mixing matrix \eqn{A} and in non-Gaussian projection pursuit this is done by first standardizaing the observations and then finding mutually orthogonal directions maximizing a convex combination of the chosen squared objective functions.

After estimation the found signals are ordered in decreasing order with respect to their objective function values.
}
\value{
A list with class 'bss' containing the following components: 
\item{W}{Estimated unmixing matrix}
\item{S}{Matrix of size \eqn{n \times k}{n x k} containing the estimated signals.}
\item{D}{Vector of the objective function values of the signals}
\item{MU}{Location vector of the data which was substracted before estimating the signal components.}
}
\references{
\cite{Virta, J., Nordhausen, K. and Oja, H., (2016), Projection Pursuit for non-Gaussian Independent Components, <https://arxiv.org/abs/1612.05445>.} 
}
\author{
Joni Virta
}



\seealso{
\code{\link{NGPPsim}, \link{NGPPest}, \link[fICA]{fICA}}
}

\examples{
# Simulated data with 2 signals

n <- 500
S <- cbind(rexp(n), runif(n),  rnorm(n))
A <- matrix(rnorm(9), ncol = 3)
X <- S \%*\% t(A)

res <- NGPP(X, 2)
res$W \%*\% A


# Iris data

X <- as.matrix(iris[, 1:4])

res <- NGPP(X, 2, nl = c("pow3", "tanh"), alpha = 0.5)
plot(res, col = iris[, 5])
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

