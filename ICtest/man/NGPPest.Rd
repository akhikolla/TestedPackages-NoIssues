\name{NGPPest}
\alias{NGPPest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Signal Subspace Dimension Testing Using non-Gaussian Projection Pursuit 
}
\description{
Estimates the dimension of the signal subspace using NGPP to conduct sequential hypothesis testing. The test statistic is a multivariate extension of the classical Jarque-Bera statistic and the distribution of it under the null hypothesis is obtained by simulation.
}
\usage{
NGPPest(X, nl = c("skew", "pow3"), alpha = 0.8, N = 500, eps = 1e-6,
        verbose = FALSE, maxiter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Numeric matrix with n rows corresponding to the observations and p columns corresponding to the variables.}
  \item{nl}{Vector of non-linearities, a convex combination of the corresponding squared objective functions of which is then used as the projection index. The choices include \code{"skew"} (skewness), \code{"pow3"} (excess kurtosis), \code{"tanh"} (\eqn{log(cosh)}) and \code{"gauss"} (Gaussian function).}
  \item{alpha}{Vector of positive weights between 0 and 1 given to the non-linearities. The length of \code{alpha} should be either one less than the number of non-linearities in which case the missing weight is chosen so that \code{alpha} sums to one, or equal to the number of non-linearities in which case the weights are used as such. No boundary checks for the weights are done.}
  \item{N}{Number of normal samples to be used in simulating the distribution of the test statistic under the null hypothesis.}
  \item{eps}{Convergence tolerance.}
  \item{verbose}{If \code{TRUE} the numbers of iterations will be printed.}
  \item{maxiter}{Maximum number of iterations.}
}
\details{
It is assumed that the data is a random sample from the model \eqn{x = m + A s} where the latent vector \eqn{s = (s_1^T, s_2^T)^T}{s = (s_1', s_2')'} consists of \eqn{k}-dimensional non-Gaussian subvector (the signal) and \eqn{p - k}-dimensional Gaussian subvector (the noise) and the components of \eqn{s} are mutually independent. Without loss of generality we further assume that the components of \eqn{s} have zero means and unit variances.

The algorithm first estimates full \eqn{p} components from the data using deflation-based NGPP with the chosen non-linearities and weighting and then tests the null hypothesis \eqn{H_0: k_{true} \leq k}{H0: k_true <= k} for each \eqn{k = 0, \ldots , p - 1}{k = 0, ..., p -1}. The testing is based on the fact that under the null hypothesis \eqn{H_0: k_{true} \leq k}{H0: k_true <= k} the distribution of the final \eqn{p - k} components is standard multivariate normal and the significance of the test can be obtained by comparing the objective function value of the \eqn{(k + 1)}th estimated components to the same quantity estimated from \code{N} samples of size \eqn{n} from \eqn{(p - k)}-dimensional standard multivariate normal distribution.

Note that if \code{maxiter} is reached at any step of the algorithm it will use the current estimated direction and continue to the next step. 
}
\value{
A list with class 'icest' containing the following components:
\item{statistic}{Test statistic, i.e. the objective function values of all estimated component.}
\item{p.value}{Obtained vector of \eqn{p}-values.}
\item{parameter}{Number \code{N} of simulated normal samples.}
\item{method}{Character string \code{"Estimation the signal subspace dimension using NGPP"}.}
\item{data.name}{Character string giving the name of the data.}
\item{W}{Estimated unmixing matrix}
\item{S}{Matrix of size \eqn{n \times p}{n x p} containing the estimated signals.}
\item{D}{Vector of the objective function values of the signals}
\item{MU}{Location vector of the data which was substracted before estimating the signal components.}
\item{conv}{Boolean vector telling for which components the algorithm converged (\code{TRUE}) and for which not (\code{FALSE}).}
}
\references{
\cite{Virta, J., Nordhausen, K. and Oja, H., (2016), Projection Pursuit for non-Gaussian Independent Components, <https://arxiv.org/abs/1612.05445>.} 
}
\author{
Joni Virta
}



\seealso{
\code{\link{NGPP}, \link{NGPPsim}}
}

\examples{
# Iris data

X <- as.matrix(iris[, 1:4])

# The number of simulations N should be increased in practical situations
# Now we settle for N = 100

res <- NGPPest(X, N = 100)
res$statistic
res$p.value
res$conv
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

