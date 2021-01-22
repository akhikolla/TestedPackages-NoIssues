\name{YPInterimTesting-package}
\alias{YPInterimTesting-package}
\docType{package}
\title{
	Interim Monitoring Using Adaptively Weighted Log-rank Test in Clinical Trials
}
\description{
This package provides monitoring boundary for interim testing using the adaptively weighted log-rank test developed by Yang and Prentice (2010). It is shown in Yang (2018) that the adaptively weighted log-rank test is asymptotically fully efficient when the limiting censoring distributions in two groups are equal. For a given spending function, the main function of the package \code{ypinterim} provides stopping boundaries at the interim looks, along with nominal p-values. The object being returned by the function \code{\link{ypinterim}} can be formatted to a table using the function \code{\link{summary}}.
}

\details{
\tabular{ll}{
Package: \tab YPInterimTesting\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2018-10-31\cr
License: \tab GPL (>= 3)\cr
}
}
\author{
Daewoo Pak and Song Yang
}
\references{
Yang, S. (2018). Interim monitoring using the adaptively weighted log-rank test in clinical trials for survival outcomes. Statistics in Medicine. <doi:10.1002/sim.7958>

Yang, S., & Prentice, R. (2010). Improved logrank-type tests for survival data using adaptive weights. Biometrics, 66(1), 30-38.
}

\keyword{Interim analysis}

\seealso{
\code{\link{ypinterim}}
}
