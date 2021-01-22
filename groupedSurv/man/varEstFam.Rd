\name{varEstFam}
\alias{varEstFam}
\title{
	Estimate the Frailty Variance for the Frailty Model
}
\description{
  A method to estimate the frailty variance for a frailty model accounting for family structure of related individuals (e.g., trios). The input data is assumed to be organized such that records for each family occur consecutively, and that records for offspring precede those for parents. The variance matrix for the random effects is assumed to be of the form \code{var*K}, where \code{K} is a matrix of kinship coefficients between family members. The following groupings are permitted: (Individual), (Offspring, Offspring), (Offspring, Parent), (Offspring, Parent, Parent), and (Offspring, Offspring, Parent, Parent). Other family structures have not been implemented. 
}

\usage{
  varEstFam(x, fam_group, fam_role, alpha, gtime, delta, lower, upper, beta=0)  	
}

\arguments{
  \item{x}{Vector of numeric variables of interest for each sample.}
  \item{fam_group}{Vector of family IDs for each sample.}
  \item{fam_role}{Vector of indicators for the role within a family of each sample, i.e., \{"Offspring","Mother","Father"\}, or \{"o","m","f"\}.}
  \item{alpha}{Vector of baseline survival rates for each time interval.}
  \item{gtime}{Vector of observed survival times for each sample.}
  \item{delta}{Vector of event indicators for each sample: 1 indicates observed event, 0 indicates censored.}
  \item{lower}{Scalar for the lower bound of the variance estimation search region.}
  \item{upper}{Scalar for the upper bound of the variance estimation search region.}
  \item{beta}{Scalar for the fixed effect size. Default is 0.}
}


\value{
  Scalar estimate of the frailty variance.
}


\examples{
# Generate dummy data	
x     <- c(0, 1, 1, 1, 2, 2, 0, 0, 0)
fam_group <- c('1', '1', '1', '2', '2', '2', '3', '3', '3')
fam_role  <- c("o","f","m","o","f","m","o","f","m")
alpha <- c(0.7500000, 0.6666667, 0.5000000, 0.0000000)
gtime <- c(1, 3, 3, 2, 1, 1, 2, 3, 1)
delta <- c(1, 0, 1, 1, 1, 0, 1, 0, 1)
lower <- 0
upper <- 2
beta  <- 1

#res   <- varEstFam(x, fam_group, fam_role, alpha, gtime, delta, lower, upper, beta) 
#res
}

