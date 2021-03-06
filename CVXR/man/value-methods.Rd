% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{value-methods}
\alias{value-methods}
\alias{value}
\alias{value<-}
\title{Get or Set Value}
\usage{
value(object)

value(object) <- value
}
\arguments{
\item{object}{A \linkS4class{Variable}, \linkS4class{Parameter}, \linkS4class{Expression}, or \linkS4class{Problem} object.}

\item{value}{A numeric scalar, vector, or matrix to assign to the object.}
}
\value{
The numeric value of the variable, parameter, or expression. If any part of the mathematical object is unknown, return \code{NA}.
}
\description{
Get or set the value of a variable, parameter, expression, or problem.
}
\examples{
lambda <- Parameter()
value(lambda)

value(lambda) <- 5
value(lambda)
}
