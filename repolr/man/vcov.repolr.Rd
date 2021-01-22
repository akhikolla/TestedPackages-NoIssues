\name{vcov.repolr}
\alias{vcov.repolr}
\title{
Calculates the Variance-Covariance Matrix for Fitted repolr Model 
}
\description{
Returns the variance-covariance matrix of the main parameters of a fitted \code{repolr} model object.
}
\usage{
\method{vcov}{repolr}(object, robust.var = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{is a model fitted using \code{repolr}.}
  \item{robust.var}{logical; if \code{TRUE}, outputs the robust variance matrix.}
  \item{\dots}{further arguments passed to or from other methods.}  
}
\details{
Default is to output the estimated robust variance matrix. However, if \code{robust.var} is set to \code{FALSE}, the naive variance matrix is reported.
}
\value{
A matrix of the estimated covariances between the parameter estimates.
}
\examples{
data(HHSpain)
mod.0 <- repolr(HHSpain~Sex*Time, data=HHSpain, categories=4, subjects="Patient",
            times=c(1,2,5), corr.mod="uniform", alpha=0.5)
vcov(mod.0, robust.var = FALSE)
}