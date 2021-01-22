\name{ virtual}
\alias{ virtual}
\docType{data}
\title{virtual Data Set for the example}
\description{
This virtual data set is created to show how to use the package.
}
\usage{data(virtual)}
\format{
  The data has the following information:
  \describe{
    \item{\code{time}}{the numeric matrix of event times for all interim looks up to the current look.}
    \item{\code{event}}{the numeric matrix of right-censoring indicators corresponding to \code{time} (event = 1, censored = 0).}
    \item{\code{group}}{the numeric vector of the group indicator (treatment = 1, control = 0).}
  }
}

\examples{
library(YPInterimTesting)
data(virtual)
time <- virtual$time
colnames(time) # Check the order of the column names.
head(time)
}
\keyword{Data}
\seealso{\code{\link{ypinterim}}}

