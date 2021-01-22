\name{polycuts}
\alias{polycuts}
\title{
Estimates Cut-point Parameters for Fitted repolr Model
}
\description{
After fitting a model using \code{repolr}, function \code{polycuts} gives estimates and standard errors for the \emph{K}-1 cut-point parameters, based on the polynomial model from the fit of \code{repolr}. Polynomial cut-point parameter estimates from the orginal model are also shown.
}
\usage{
polycuts(object, digits = 3, robust.var = TRUE)
}
\arguments{
  \item{object}{is a model fitted using \code{repolr}.}
  \item{digits}{the number of decimal places to display in reported summaries.}
  \item{robust.var}{a logical variable: if \code{TRUE} standard errors are based on robust variance estimates, otherwise naive estimates are used.}  
}
\value{
 \item{coef}{polynomial parameter estimates from \code{repolr}.}
 \item{poly}{a vector of \emph{K}-1 cut-point parameters.}
 \item{order}{the order of the polynomial.} 
}
\keyword{models}

\examples{

data(HHSpain)
mod.0 <- repolr(HHSpain~Sex*Time, data=HHSpain, categories=4, subjects="Patient",
            times=c(1,2,5), corr.mod="uniform", alpha=0.5)
summary(mod.0)
mod.1 <- update(mod.0, poly=1)
summary(mod.1)
polycuts(mod.1)
mod.2 <- update(mod.0, poly=2)
summary(mod.2)
polycuts(mod.2)

}