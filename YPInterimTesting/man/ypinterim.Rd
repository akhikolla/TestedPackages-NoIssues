\name{ypinterim}
\alias{ypinterim}
\alias{ypinterim.default}
\alias{summary.ypinterim}
\alias{print.ypinterim}
\title{The main function of the package performing the adaptively weighted log-rank test for interim monitoring}
\usage{
\method{ypinterim}{default}(time, event, group, spendfun, critvalue = NULL,
  repnum = 1E4, bound = 50, seed.fix = 0, ...)
}
\arguments{
\item{...}{for S4 method only.}

\item{time}{A numeric matrix of event times for all interim looks up to the current look. The columns of the matrix represent interim looks in date order.}

\item{event}{A numeric matrix of right-censoring indicators corresponding to \code{time} (event = 1, censored = 0). The dimension of the matrix should be equivalent to that of \code{time}.}

\item{group}{A numeric vector of the group indicator (treatment = 1, control = 0).}

\item{spendfun}{A numeric vector of the values for the spending function \eqn{\alpha(t)}. Must be user-specified.}

\item{critvalue}{A vector of the previous critical values for the prior interim looks. If it is not supplied, the critical values at all the interim looks will be obtained all at once using the re-sampling method. This is only for the scenario when the user wants to test with historical data where all interim data are available. When using the package for an ongoing trial, at the first look, when the user supplies the data and the value of the spending function at the first look, the package provides the stopping boundary, the test statistic value and the nominal p-values at the first look. Then, at the \eqn{i}-th look, \eqn{i>1}, when the user supplies the data and the values of the spending function at the \eqn{i} looks, as well as the previous \eqn{i-1} critical values that have already been obtained, the package provides the stopping boundary, the test statistics value and the nominal p-values at the \eqn{i}th look, and display results for all \eqn{i} looks.}

\item{repnum}{The number of replications for the re-sampling method used in the calculations. The default value is one million. When the  alpha allocated is less than 1E-4 at some early looks, the package uses numerical integration  to obtain the critical values at those looks.}
  
\item{bound}{
A boundary (\code{-bound}, \code{bound}) for estimating the parameters in the short-term and long-term hazard ratio model (Yang and Prentice, 2005). These parameters are \eqn{\beta_1} and \eqn{\beta_2} in their notations. The default boundary is \eqn{(-50, 50)}.}

\item{seed.fix}{A value for \code{set.seed}. The default value is 0. Should be FALSE or a numeric number.}
}

\description{
Provides stopping boundaries and nominal p-values at the interim looks using the adaptively weighted log-rank test developed by Yang and Prentice (2010). The full details are available on Yang (2018).
}
\details{
The object being returned by the function \code{\link{ypinterim}} can be formatted to a table using the function \code{\link{summary}}.
}

\references{
Yang, S. (2018). Interim monitoring using the adaptively weighted log-rank test in clinical trials for survival outcomes. Statistics in Medicine. <doi:10.1002/sim.7958>

Yang, S., & Prentice, R. (2010). Improved logrank-type tests for survival data using adaptive weights. Biometrics, 66(1), 30-38.

Yang, S., & Prentice, R. (2005). Semiparametric analysis of short-term and long-term hazard ratios with two-sample survival data. Biometrika, 92(1), 1-17.
}

\examples{
library(YPInterimTesting)
data(virtual)

time <- virtual$time
event <- virtual$event
group <- virtual$group

spendfun <- c(1.3E-5, 4.4E-4, 0.003, 0.008)

result_all <- ypinterim(time, event, group, spendfun=spendfun)
result_all
summary(result_all)

## The example shows how to test the package with a historical data where 
## interim data at all looks are available.

## When calculating the boundary at the first look,
## the spending function value at the first look is needed:

time <- virtual$time[, 1]
event <- virtual$event[, 1]
group <- virtual$group

spendfun <- c(1.3E-5)

result_look1 <- ypinterim(time, event, group, spendfun=spendfun)
result_look1
summary(result_look1)

## When calculating the boundary at the second look,
## the spending function at the two looks,
## and boundary at the first look, should be supplied:

time <- virtual$time[, 1:2]
event <- virtual$event[, 1:2]
group <- virtual$group

spendfun <- c(1.3E-5, 4.4E-4)
critvalue <- c(4.36) # the boundary of the first look is supplied.

result_look2 <- ypinterim(time, event, group, spendfun=spendfun, critvalue = critvalue)
result_look2
summary(result_look2)

## Similarly, when calculating the boundary at the third look,
## the spending function at the three looks, and boundaries
## at the first two looks, should be supplied:
time <- virtual$time[, 1:3]
event <- virtual$event[, 1:3]
group <- virtual$group

spendfun <- c(1.3E-5, 4.4E-4, 0.003)
critvalue <- c(4.36, 3.42) # the boundaries at the first two looks are supplied.

result_look3 <- ypinterim(time, event, group, spendfun=spendfun, critvalue = critvalue)
result_look3
summary(result_look3)

}

\keyword{ypinterim}