\name{densEstBayes.control}
\alias{densEstBayes.control}
\title{Controlling density estimation via Bayesian inference engines}
\description{Function for optional use in calls to \code{densEstBayes()} to control convergence values and other specifications for density estimation using Bayesian inference engines.}
\usage{
densEstBayes.control(range.x = NULL,numBins = 401,numBasis = 50,sigmabeta = 1e5,
                     ssigma = 1000,convToler = 1e-5,maxIter = 500,nWarm = NULL,
                     nKept = NULL,nThin = 1,msgCode = 1)
}
\arguments{
\item{range.x}{A two-component array such that the density estimate is obtained over the interval between range.x[1] and range.x[2]. The default value of range.x[1] is 1.05*min(x)-0.05*max(x) and the default value of range.x[2] is 1.05*max(x)-0.05*min(x).}
\item{numBins}{The number of bins used for binning of the data. The default is 401.}
\item{numBasis}{The number of cubic O'Sullivan spline basis functions. The default is 50.}
\item{sigmabeta}{The prior standard deviation of the fixed effects coefficients corresponding to the linear component of the fit. The default is 100000.}
\item{ssigma}{The prior scale parameter of the standard deviation parameter. The default is 1000.}
\item{convToler}{The convergence tolerance value used for determination of convergence of when the method is semiparametric mean field variational Bayes. Convergence is deemed to have occured when the relative change in the q-density expectation of the reciprocal variance parameter is below convToler. The default is 0.00001.}
\item{maxIter}{The maximum number of iterations when the method is mean field variational Bayes. The default is 500.}
\item{nWarm}{The size of the Markov chain Monte Carlo warmup, a positive integer, when the method is Markov chain Monte Carlo. The default for method = "HMC" and method = "NUTS" is 1000 and the default for method = "slice" is 100.}
\item{nKept}{The size of the kept Markov chain Monte Carlo samples, a positive integer, when the method is Markov chain Monte Carlo. The default is 1000.}
\item{nThin}{The thinning factor for the kept Markov chain Monte Carlo samples, a positive integer, when the method is Markov chain Monte Carlo. The default is 1.}
\item{msgCode}{Code for the nature of progress messages printed to the screen when the method is either "HMC", "NUTS" or "sl
ice". If msgCode=0 then no progress messages are printed. If msgCode=1 then a messages printed each time approximately each 10\% of the sampling is completed. If msgCode=2 then a messages printed each time approximately each 10\% of the sampling is completed. The default is 1.}
}

\value{A list containing values of each of the fifteen control parameters, packaged to supply the \code{control} argument to \code{densEstBayes}. The values for \code{densEstBayes.control} can be specified in the call to \code{densEstBayes}.
}
\author{Matt P. Wand \email{matt.wand@uts.edu.au}}
\references{Wand, M.P. and Yu, J.C.F. (2020). 
Density estimation via Bayesian inference engines. Submitted.}
\examples{
library(densEstBayes) ; set.seed(1)
x <- rMarronWand(1000,8)
dest <- densEstBayes(x,method = "SMFVB")
plot(dest) ; rug(x,col="dodgerblue",quiet = TRUE)
xg <- seq(-3,3,length = 1001)
trueDensg <- dMarronWand(xg,8)
lines(xg,trueDensg,col = "indianred3")

# Now modify some of the control values:

destControlled <- densEstBayes(x,method = "SMFVB",control = densEstBayes.control(numBins = 201,
                               numBasis = 35,sigmabeta = 1000,ssigma = 100,convToler = 1e-4))
plot(destControlled) ; rug(x,col = "dodgerblue",quiet = TRUE)
lines(xg,trueDensg,col = "indianred3")
}



