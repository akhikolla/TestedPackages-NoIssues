\name{betaEst}
\alias{betaEst}
\title{
  Estimate the Parameter of Interest for the Grouped Survival Model
}
\description{
  A method to estimate the parameter of interest given the baseline survival rates for each time interval and optional covariate nuisance parameters for a grouped survival model. The function supports zero, single, or multiple covariates as dictated by the input data. 
}

\usage{
  betaEst(x, Z=NULL, alpha, theta=NULL, gtime, delta)  	
}

\arguments{
  \item{x}{Vector of numeric variable of interest for each sample.}
  \item{Z}{Optional \code{data.frame} or \code{matrix} of numeric covariate values for each sample.}
  \item{alpha}{Vector of baseline survival rates for each time interval.}
  \item{theta}{Optional vector of estimated nuisance parameters for the covariates.}
  \item{gtime}{Vector of observed survival times for each sample.}
  \item{delta}{Vector of event indicators for each sample: 1 indicates observed event, 0 indicates censored.}
}


\value{
  Scalar estimate of the parameter of interest.
}


\examples{
# Generate dummy data	
x     <- c(1, 0, 1, 1, 1, 0, 1, 0, 1)
cov1  <- c(1, 2, 2, 2, 1, 1, 0, 1, 1)
cov2  <- c(2, 2, 1, 0, 1, 0, 1, 1, 0.5) 
Z     <- cbind(cov1, cov2)
gtime <- c(1, 3, 3, 2, 1, 1, 2, 3, 1)
delta <- c(1, 0, 1, 1, 1, 0, 1, 0, 1)

#theta <- thetaEst(Z, gtime, delta) 

#res   <- betaEst(x, Z, theta$alpha, theta$theta, gtime, delta) 
#res
}

