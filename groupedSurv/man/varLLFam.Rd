\name{varLLFam}
\alias{varLLFam}
\title{
	Log-Likelihood of the Frailty Variance for the Discrete-Time Frailty Model
}
\description{
  A method to compute the log-likelihood of the frailty variance for a frailty model accounting for family structure of related individuals (e.g., trios). The input data is assumed to be organized such that records for each family occur consecutively, and that records for offspring precede those for parents. The variance matrix for the random effects is assumed to be of the form \code{var*K}, where \code{K} is a matrix of kinship coefficients between family members. The following groupings are permitted: (Individual), (Offspring, Offspring), (Offspring, Parent), (Offspring, Parent, Parent), and (Offspring, Offspring, Parent, Parent). Other family structures have not been implemented. 
}

\usage{
  varLLFam(x, fam_group, fam_role, alpha, var, gtime, delta, beta=0)  	
}

\arguments{
  \item{x}{Vector of numeric variables of interest for each sample.}
  \item{fam_group}{Vector of family IDs for each sample.}
  \item{fam_role}{Vector of indicators for the role within a family of each sample, i.e., \{"Offspring","Mother","Father"\}, or \{"o","m","f"\}.}
	\item{alpha}{Vector of baseline survival rates for each time interval.}
  \item{var}{Scalar for frailty variance.}
  \item{gtime}{Vector of observed survival times for each sample.}
  \item{delta}{Vector of event indicators for each sample: 1 indicates observed event, 0 indicates censored.}
  \item{beta}{Scalar for the fixed effect size. Default is 0.}
}

\value{
   Scalar for the log-likelihood of the frailty variance.
}


\examples{
# Generate dummy data	
x     <- c(0, 1, 1, 1, 2, 2, 0, 0, 0)
fam_group <- c('1', '1', '1', '2', '2', '2', '3', '3', '3')
fam_role <- c("o","f","m","o","f","m","o","f","m")
alpha <- c(0.7500000, 0.6666667, 0.5000000, 0.0000000)
var <- 0.1
gtime <- c(1, 3, 3, 2, 1, 1, 2, 3, 1)
delta <- c(1, 0, 1, 1, 1, 0, 1, 0, 1)
beta  <- 1.0

#res   <- varLLFam(x, fam_group, fam_role, alpha, var, gtime, delta, beta) 
#res
}

