\name{summary.repolr}
\alias{summary.repolr}
\alias{print.summary.repolr}
\title{
Summary of Fitted repolr Model
}
\description{
Function to summarise the fit of a repolr model.
}
\usage{
\method{summary}{repolr}(object, digits, robust.var, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{fitted model.}
  \item{digits}{integer for number formatting.}
  \item{robust.var}{logical; if \code{TRUE}, standard errors are based on robust variance estimates}
  \item{\dots}{further arguments passed to or from other methods.}  
}
\details{
Default is to use robust variance estimates. However, if \code{robust.var} is set to \code{FALSE}, naive variance estimates are used.
}
