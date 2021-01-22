\name{mobility}
\alias{mobility}
\docType{data}
\title{Mobility after Hip Fracture Fixation Surgery}
\description{
Patient mobility score after hip fracture fixation surgery, using two distinct procedures. 
}
\usage{data(mobility)}
\format{
A data frame with 600 observations on the following variables.
  \describe{
    \item{\code{subject}}{patient identifier.}
    \item{\code{time}}{assessemnt occasions; month 1, 2, 3 and 4.}
    \item{\code{treat}}{intervention group; \code{A} or \code{B}.}
    \item{\code{age}}{patient age in years.}
    \item{\code{gender}}{gender; \code{F} or \code{M}.}
    \item{\code{mobility}}{mobility score; 1, 2 and 3.}
  }
}
\examples{
data(mobility)
head(mobility)
}
\keyword{datasets}
