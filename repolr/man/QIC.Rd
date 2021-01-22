\name{QIC}
\alias{QIC}
\title{Quasilikelihood Information Criterion}
\description{
The quasilikelihood information criterion (QIC) developed by Pan (2001) is a modification of the Akaike information criterion (AIC) for models fitted by GEE. \code{QIC} is used for choosing the best correaltion structure and \code{QICu} is used for choosing the best subset of covariates. The quasilikelihood (\code{QLike}) is also reported for completeness. When choosing between two or more models, with different subset of covariates, the one with the smallest \code{QICu} measure is preferred and similarly, when choosing between competing correlation structures, with the same subset of covariates in both, the model with the smallest \code{QIC} measure is preferred.
}
\usage{
QIC(object, digits = 3)
}
\arguments{
  \item{object}{is a fitted model using \code{repolr}.}
  \item{digits}{the number of decimal places to display in reported summaries.}
}
\value{
 \item{QLike}{model quasilikelihood.}
 \item{QIC}{model \code{QIC}.}
 \item{QICu}{model \code{QICu}.} 
}
\references{
Pan W. Akaikes information criterion in generalized estimating equations. \emph{Biometrics} 2001; 57:120-125.\cr
}
\keyword{models}

\examples{

data(HHSpain)
mod.0 <- repolr(HHSpain~Time, data=HHSpain, categories=4, subjects="Patient",
            times=c(1,2,5), corr.mod="independence", alpha=0.5)
QIC(mod.0)
QIC(update(mod.0, formula = HHSpain~Time + Sex))$QICu
QIC(update(mod.0, formula = HHSpain~Time * Sex))$QICu

}

