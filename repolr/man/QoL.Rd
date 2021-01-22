\name{QoL}
\alias{QoL}
\docType{data}
\title{Quality of Life Scores}
\description{
Quality of life scores after hip replacement.
}
\usage{data(QoL)}
\format{
  A data frame with 336 observations on the following four variables.
  \describe{
    \item{\code{QoL}}{a numeric vector.}
    \item{\code{Patient}}{a numeric vector.}
    \item{\code{Time}}{a numeric vector.}
    \item{\code{Treat}}{a factor with levels \code{A} \code{B}.}
  }
}

\examples{
data(QoL)
head(QoL)
}
\keyword{datasets}
