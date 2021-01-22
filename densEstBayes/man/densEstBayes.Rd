\name{densEstBayes}
\alias{densEstBayes}
\title{Density estimation via Bayesian inference engines}
\description{Construction of a high quality density estimate from a random sample is a fundamental problem in Statistics. This package delivers solutions to this problem by calling upon Bayesian inference engine algorithms with roots in the Machine Learning literature. Both stochastic and faster deterministic algorithms are provided. The stochastic options are Hamiltonian Monte Carlo and the no U-turn sampler (through the Stan inference engine) and slice sampling. The deterministic option is semiparametric mean field variational Bayes. The last of these options provides a very quick Bayesian density estimate, whereas the other options take longer to compute but, according to extensive simulation testing, are highly accurate. The essence of the approach is conversion of the density estimation problem to a Bayesian Poisson mixed model problem via binning on a fine grid. An attractive by-product of the binning approach is the ability to handle very large sample sizes. The number of bins, rather than the sample size, is the limiting factor.}
\usage{
densEstBayes(x,method,verbose,control)
}
\arguments{
\item{x}{Vector containing a continuous univariate data set.}
\item{method}{Character string for specifying the method to be used:\cr
  "HMC" = Hamiltonian Monte Carlo,\cr
  "NUTS" = the no U-turn sampler,\cr
  "slice" = slice sampling,\cr
  "SMFVB" = semiparametric mean field variational Bayes.}
\item{verbose}{Boolean variable for specifying whether or not progress messages are printed to the console. The default is TRUE.}
\item{control}{Function for controlling Markov chain Monte Carlo sample sizes and other specifications.}
}
\details{The crux of the Bayesian density estimation approach is to bin the input data on a fine equally-spaced grid and treat the bin counts as response data in a Bayesian Poisson mixed model-based penalized splines model. Section 8 of Eilers and Marx (1996) provides details on couching the density estimation problem as a Poisson penalized splines problem.\cr

Fitting and inference for the resultant Bayesian Poisson mixed model is carried out using one of four possible approaches as specified by the \code{method} argument; one of which is deterministic ("SMFVB") and three of which are stochastic ("HMC", "NUTS" and "slice").

The settings \code{method} = "HMC" and \code{method} = "NUTS" correspond, respectively, to the Markov chain Monte Carlo approaches Hamiltonian Monte Carlo (e.g. Neal, 2011) and the no-U-turn sampler (Hoffman and Gelman, 2014) and fitting and inference is performed using the Stan Bayesian inference engine via the \code{rstan} package (Stan Development Team, 2017). 

The setting \code{method} = "slice" corresponds to slice sampling (e.g. Neal, 2003).

The setting \code{method} = "SMFVB" uses a semiparametric mean field variational Bayes approach as summarised in Algorithm 1 of Luts and Wand (2015). Comparatively speaking, this method leads to a very quick density estimate but at some accuracy and convergence reliability costs.
}
\value{An object of class \code{densEstBayes}, which is a list with the following components:
\item{method}{the value of \code{method}} 
\item{range.x}{a numerical vector with 2 entries such that range.x[1] is the lower limit and range.x[2] is the upper limit of the interval over which the density estimate was computed.}
\item{intKnots}{a numerical vector containing the locations of the interior knots used in the cubic O'Sullivan spline basis}
\item{determFitObj}{if \code{method} is "SMFVB" then \code{determFitObj} is a list with the following components:\cr
mu.q.betau: numerical vector containing the mean vector of the Multivariate Normal approximate posterior density function of the linear component coefficients (first two entries) and penalised spline coefficients (remainining entries).\cr
Sigma.q.betau: numerical matrix containing the covariance matrix of the Multivariate Normal approximate posterior density function of linear and penalised spline coefficients vector.\cr 
kappa.q.sigsq: numerical scalar containing the shape parameter of the Inverse Gamma approximate posterior density function of the variance parameter in the mixed model-based penalised spline model.\cr
lambda.q.sigsq: numerical scalar containing the scale parameter of the Inverse Gamma approximate posterior density function of the variance parameter in the mixed model-based penalised spline model.\cr
If \code{method} is "HMC", "NUTS" or "slice" then \code{determFitObj} is returned as NULL.
}
\item{stochaFitObj}{if \code{method} is "HMC", "NUTS" or "slice" then \code{stochaFitObj} is a list with the following components:\cr
betauMCMC: numerical matrix containing Markov chain Monte Carlo draws from the posterior distribution of the basis function coefficients, with rows corresponding to the linear component coefficients (first two rows) and penalised spline coefficients (remainining rows) and columns corresponding to the samples.\cr
sigsqMCMC: numerical vector containing Markov chain Monte Carlo draws from the posterior distribution of the variance parameter in the mixed model-based penalised spline model.\cr
If \code{method} is "SMFVB" then \code{stochaFitObj} is returned as NULL.}
\item{xname}{ONE OF THE LAST THINGS TO ADD.}
\item{sampHexTran}{Sample hexiles of the transformed input data to be used for possible checking Markov chains Monte Carlo samples.}
}
\author{Matt P. Wand \email{matt.wand@uts.edu.au}}
\references{
Eilers, P.H.C. and Marx, B.D. (1996).
Flexible smoothing with B-splines and penalties (with discussion).
\emph{Statistical Science}, \bold{11}, 89-121.

Hoffman, M.D. and Gelman, A. (2014). The no-U-turn sampler: adaptively setting path lengths in Hamiltonian Monte Carlo.
\emph{Journal of Machine Learning Research}, \bold{15}, 1593-1623.

Minka, T. (2005). Divergence measures and message passing. \emph{Microsoft Research Technical Report Series}, 
\bold{MSR-TR-2005-173}, 1-17.

Luts, J. and Wand, M.P. (2015). Variational inference for count response semiparametric regression.
\emph{Bayesian Analysis}, \bold{10}, 991-1023.

Neal, R. (2003). Slice sampling (with discussion). \emph{The Annals of Statistics}, \bold{31}, 705-767.

Neal, R. (2011). MCMC using Hamiltonian dynamics. In \emph{Handbook of Markov Chain Monte Carlo}, eds. S. Brooks, A. Gelman, G.L. Jones and X.-L. Meng. Boca Raton, Florida: Chapman & Hall/CRC Press.

Stan Development Team. (2020). \emph{Stan Modeling Language
User's Guide and Reference Manual}, \url{https://mc-stan.org}.

Wainwright, M.J. and Jordan, M.I. (2008). Graphical models, exponential families, and variational inference.
\emph{Foundations and Trends in Machine Learning}, \bold{1}, 1-305.

Wand, M.P. and Yu, J.C.F. (2020). Density estimation via Bayesian inference engines. Submitted.
}
\examples{
library(densEstBayes) ; set.seed(1)
x <- rMarronWand(1000,8)
dest <- densEstBayes(x,method = "SMFVB")
plot(dest) ; rug(x,col = "dodgerblue",quiet = TRUE)
xg <- seq(-3,3,length = 1001)
trueDensg <- dMarronWand(xg,8)
lines(xg,trueDensg,col = "indianred3")
\donttest{
if (require("rstan"))
{
   dest <- densEstBayes(x,method = "NUTS")
   plot(dest) ; rug(x,col = "dodgerblue")
   xg <- seq(-3,3,length = 1001)
   trueDensg <- dMarronWand(xg,8)
   lines(xg,trueDensg,col = "indianred3")
}
}
}


