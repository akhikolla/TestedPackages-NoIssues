\name{confint.repolr}
\alias{confint.repolr}
\title{
Confidence Intervals for repolr Model Parameters
}
\description{
Computes confidence intervals for one or more parameters in a fitted \code{repolr} model object.
}
\usage{
\method{confint}{repolr}(object, parm, level = 0.95, robust.var = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{is a model fitted using \code{repolr}.}
  \item{parm}{a specification of which parameters are to be used, either a vector of numbers or a vector of names. If missing, all parameters are considered.}
  \item{level}{the confidence level required.}
  \item{robust.var}{logical; if \code{TRUE}, intervals are based on the robust variance matrix.}
  \item{\dots}{further arguments passed to or from other methods.}  
}
\details{
The method assumes normality and uses as default the estimated robust variance matrix.
}
\value{
A matrix (or vector) with columns giving lower and upper confidence limits for each parameter.
}
\examples{
data(HHSpain)
mod.0 <- repolr(HHSpain~Sex*Time, data=HHSpain, categories=4, subjects="Patient",
            times=c(1,2,5), corr.mod="uniform", alpha=0.5)
confint(mod.0, robust.var = FALSE)
}