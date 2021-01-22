\name{work.corr}
\alias{work.corr}
\title{
Working Correlation Matrix for Fitted repolr Model
}
\description{
The working correlation matrix for the selected model; \dQuote{\code{ar1}}, \dQuote{\code{uniform}} or \dQuote{\code{independence}}.
}
\usage{
work.corr(object, digits = 3)
}
\arguments{
  \item{object}{is a model fitted using \code{repolr}.}
  \item{digits}{integer for number formatting.}
}

\value{
 A \emph{T}(\emph{K}-1) correlation matrix.
}

\examples{
data(HHSpain)
mod.0 <- repolr(HHSpain~Sex*Time, data=HHSpain, categories=4, subjects="Patient",
            times=c(1,2,5), corr.mod="uniform", alpha=0.5)
work.corr(mod.0)
}