\name{thetaEst}
\alias{thetaEst}
\title{
  Estimate the Baseline Survival Rates and Covariate Parameters for the Grouped Survival Model
}
\description{
  A method to estimate the baseline survival rate for each time interval and the covariate nuisance parameters for a grouped survival model. The estimation is conducted under the null hypothesis, i.e., that there is no effect from the variable of interest.
}

\usage{
  thetaEst(Z=NULL, gtime, delta, method="BFGS")  	
}

\arguments{
  \item{Z}{Optional \code{data.frame} or \code{matrix} of numeric covariate values for each sample.}
  \item{gtime}{Vector of observed survival times for each sample.}
  \item{delta}{Vector of event indicators for each sample: 1 indicates observed event, 0 indicates censored.}
  \item{method}{String indicating optimization method, passed to \code{stats::optim}. Supports \code{"BFGS"} and \code{"CG"}. Default is \code{"BFGS"}.}
}


\value{
  A list containing two vectors: \code{alpha}, the baseline survival rate for each time interval, and \code{theta}, the estimated nuisance parameters of the covariates.
}


\examples{
# Generate dummy data	
cov1  <- c(1, 2, 2, 2, 1, 1, 0, 1, 1)
cov2  <- c(2, 2, 1, 0, 1, 0, 1, 1, 0.5) 
Z     <- cbind(cov1, cov2)
gtime <- c(1, 3, 3, 2, 1, 1, 2, 3, 1)
delta <- c(1, 0, 1, 1, 1, 0, 1, 0, 1)

#res   <- thetaEst(Z, gtime, delta) 
#res
}

