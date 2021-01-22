\name{VARfit-methods}
\alias{coef.VARfit}
\alias{print.VARfit}
\alias{residuals.VARfit}

\title{
Methods for class VARfit
}
\description{
\code{residuals.VARfit()} returns the residuals and \code{coef.VARfit()} returns the coefficients of a fitted VAR model of class 'VARfit', while \code{print.VARfit()} prints the parameters of the estimated model.
}

\usage{
\method{coef}{VARfit}(object, \dots)
\method{print}{VARfit}(x, \dots)
\method{residuals}{VARfit}(object, \dots)
}

\arguments{
  \item{x, object}{
the fitted VAR model of class 'VARfit' (as returned by the function \code{\link{VARfit}}).
}
  \item{\dots}{
  additional arguments.
}
}