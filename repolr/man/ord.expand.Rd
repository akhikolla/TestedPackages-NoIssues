\name{ord.expand}
\alias{ord.expand}
\title{Expand Ordinal Data for Model Fitting with repolr}
\description{
Expands ordinal score data into an appropriate form, \emph{K}-1 new binary variables for an ordinal score with \emph{K} categories, for model fitting using \code{repolr}. 
}
\usage{
ord.expand(space, formula, times, poly, data, subjects, categories)
}
\arguments{
  \item{space}{a vector indicating the category spacing when fitting the polynomial model.}
  \item{formula}{a formula, as for other regression models.}
  \item{times}{a vector of times which occur within subject clusters.}
  \item{poly}{a numeric variable indicating the order of the polynomial contrasts used for the cut-point model.}  
  \item{data}{a dataframe in which to interpret the variables occurring in the formula.}
  \item{subjects}{a character string specifying the name of the subject variable.}
  \item{categories}{a numeric variable indicating the number of ordinal score categories.}
}
\details{
For internal use with \code{repolr}.
}
\value{
An expanded dataframe.
}
\keyword{models}
