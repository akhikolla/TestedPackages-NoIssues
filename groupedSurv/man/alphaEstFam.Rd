\name{alphaEstFam}
\alias{alphaEstFam}
\title{
	Estimate the Interval Baseline Survival Rates for the Grouped Failure Time Model
}
\description{
  A method to estimate the baseline survival rate for each time interval for a grouped failure time model. The estimation is conducted under the null hypothesis, i.e., that there is no effect from the variable of interest, and is naive to family structure.
}

\usage{
  alphaEstFam(gtime, delta)  	
}

\arguments{
  \item{gtime}{Vector of observed survival times for each sample.}
  \item{delta}{Vector of event indicators for each sample: 1 indicates observed event, 0 indicates censored.}
}


\value{
  Vector of estimates of the baseline survival rates for each time interval.
}


\examples{
# Generate dummy data	
gtime <- c(1, 3, 3, 2, 1, 1, 2, 3, 1)
delta <- c(1, 0, 1, 1, 1, 0, 1, 0, 1)

#res   <- alphaEstFam(gtime, delta) 
#res
}

