\name{NGPPsim}
\alias{NGPPsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Signal Subspace Dimension Testing Using non-Gaussian Projection Pursuit 
}
\description{
Tests whether the true dimension of the signal subspace is less than or equal to a given \eqn{k}. The test statistic is a multivariate extension of the classical Jarque-Bera statistic and the distribution of it under the null hypothesis is obtained by simulation.
}
\usage{
NGPPsim(X, k, nl = c("skew", "pow3"), alpha = 0.8, N = 1000, eps = 1e-6,
        verbose = FALSE, maxiter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Numeric matrix with n rows corresponding to the observations and p columns corresponding to the variables.}
  \item{k}{Number of components to estimate, \eqn{1 \leq k \leq p}{1 <= k <= p}.}
  \item{nl}{Vector of non-linearities, a convex combination of the corresponding squared objective functions of which is then used as the projection index. The choices include \code{"skew"} (skewness), \code{"pow3"} (excess kurtosis), \code{"tanh"} (\eqn{log(cosh)}) and \code{"gauss"} (Gaussian function).}
  \item{alpha}{Vector of positive weights between 0 and 1 given to the non-linearities. The length of \code{alpha} should be either one less than the number of non-linearities in which case the missing weight is chosen so that \code{alpha} sums to one, or equal to the number of non-linearities in which case the weights are used as such. No boundary checks for the weights are done.}
  \item{N}{Number of normal samples to be used in simulating the distribution of the test statistic under the null hypothesis.}
  \item{eps}{Convergence tolerance.}
  \item{verbose}{If \code{TRUE} the numbers of iterations will be printed.}
  \item{maxiter}{Maximum number of iterations.}
}
\details{
It is assumed that the data is a random sample from the model \eqn{x = m + A s} where the latent vector \eqn{s = (s_1^T, s_2^T)^T}{s = (s_1', s_2')'} consists of \eqn{k}-dimensional non-Gaussian subvector (the signal) and \eqn{p - k}-dimensional Gaussian subvector (the noise) and the components of \eqn{s} are mutually independent. Without loss of generality we further assume that the components of \eqn{s} have zero means and unit variances.

To test the null hypothesis \eqn{H_0: k_{true} \leq k}{H0: k_true <= k} the algorithm first estimates \eqn{k + 1} components using delfation-based NGPP with the chosen non-linearities and weighting. Under the null hypothesis the distribution of the final \eqn{p - k} components is standard multivariate normal and the significance of the test is obtained by comparing the objective function value of the \eqn{(k + 1)}th estimated components to the same quantity estimated from \code{N} samples of size \eqn{n} from \eqn{(p - k)}-dimensional standard multivariate normal distribution.

Note that if \code{maxiter} is reached at any step of the algorithm it will use the current estimated direction and continue to the next step. 
}
\value{
A list with class 'ictest', inheriting from the class 'hctest', containing the following components:
\item{statistic}{Test statistic, i.e. the objective function value of the (\code{k} + 1)th estimated component.}
\item{p.value}{Obtained \eqn{p}-value.}
\item{parameter}{Number \code{N} of simulated normal samples.}
\item{method}{Character string denoting which test was performed.}
\item{data.name}{Character string giving the name of the data.}
\item{alternative}{Alternative hypothesis, i.e. \code{"There are less than p - k Gaussian components"}.}
\item{k}{Tested dimension \code{k}.}
\item{W}{Estimated unmixing matrix}
\item{S}{Matrix of size \eqn{n \times (k + 1)}{n x (k + 1)} containing the estimated signals.}
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
\code{\link{NGPP}, \link{NGPPest}}
}

\examples{
# Simulated data with 2 signals and 2 noise components

n <- 200
S <- cbind(rexp(n), rbeta(n, 1, 2), rnorm(n), rnorm(n))
A <- matrix(rnorm(16), ncol = 4)
X <- S \%*\% t(A)

# The number of simulations N should be increased in practical situations
# Now we settle for N = 100

res1 <- NGPPsim(X, 1, N = 100)
res1
screeplot(res1)

res2 <- NGPPsim(X, 2, N = 100)
res2
screeplot(res2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

