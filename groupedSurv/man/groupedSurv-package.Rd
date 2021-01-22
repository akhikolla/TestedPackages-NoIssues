\name{groupedSurv-package}
\alias{groupedSurv-package}
\docType{package}
\title{
	Efficient Estimation of Grouped Survival Models
}
\description{
These \code{Rcpp}-based functions compute the efficient score statistics for grouped time-to-event data (Prentice and Gloeckler, 1978), with the optional inclusion of baseline covariates. Functions for estimating the parameter of interest and nuisance parameters, including baseline hazards, using maximum likelihood are also provided. A parallel set of functions allow for the incorporation of family structure of related individuals (e.g., trios). Note that the current implementation of the frailty model (Ripatti and Palmgren, 2000) is sensitive to departures from model assumptions, and should be considered experimental. For these data, the exact proportional-hazards-model-based likelihood is computed by evaluating multiple variable integration. The integration is accomplished using the \code{Cuhre} algorithm from the \code{Cuba} library (Hahn, 2005), and the source files of the \code{Cuhre} function are included in this package. The maximization process is carried out using Brent's algorithm, with the \code{C++} code file from John Burkardt and John Denker (Brent, 2002). 

License: GPL (>= 2)}
\details{
\tabular{ll}{
Package: \tab groupedSurv\cr
Type: \tab Package\cr
Version: \tab 1.0.4.2\cr
Date: \tab 2020-11-11\cr
License: \tab GPL-3\cr
}
Please refer to the individual function documentation or the included vignette for more information. The package vignette serves as a tutorial for using this package.
}

\author{
  Jiaxing Lin \email{jiaxing.lin@duke.edu}, Alexander Sibley, Tracy Truong, Kouros Owzar, Zhiguo Li; Contributors: Yu Jiang, Janice McCarthy, Andrew Allen
}

\references{
  Prentice, R.L. and Gloeckler, L.A. (1978). Regression analysis of grouped survival data with application to breast cancer data. Biometrics, 34:1, 57-67.\cr
	Ripatti, S. and Palmgren, J. (2000). Estimation of Multivariate Frailty Models Using Penalized Partial Likelihood. Biometrics, 56, 1016-1022.\cr
	Hahn, T. (2005). Cuba-a library for multidimensional numerical integration, Computer Physics Communications, 168, 78-95.\cr
	Brent, R. (2002). Algorithms for Minimization without Derivatives. Dover, ISBN 0-486-41998-3
}

\keyword{ package }



