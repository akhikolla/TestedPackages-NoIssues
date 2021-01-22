\name{eBsc-package}
\alias{eBsc-package}
\docType{package}
\title{
Empirical Bayes Smoothing Splines with Correlated Errors
}
\description{
Empirical Bayes smoothing splines with correlated errors. The method uses a recursive algorithm for signal extraction with a non-parametric estimation of the correlation matrix of the errors.
}
\details{
\tabular{ll}{
Package: \tab eBsc\cr
Version: \tab 4.11\cr
Date: \tab 2020-08-21\cr
Depends: \tab Brobdingnag, parallel, nlme, Matrix, MASS\cr
}

Index:
\preformatted{
eBsc                  Empirical Bayes smoothing splines with correlated errors
plot.eBsc             Plots fitted curves from the filter
summary.eBsc          Summary information of the error
}
The function \code{eBsc()} is used to fit the model. Using the resulting 
\code{eBsc} object and summary information on the errors can be printed using \code{\link{summary}}.
}

\author{
Francisco Rosales, Paulo Serra, Tatyana Krivobokova
Maintainer: Francisco Rosales <francisco.rosales-marticorena@protonmail.com>
}
\examples{
# simulated data for non-correlated errors
library(eBsc)
n <- 250
sigma <- 0.05
beta <- function(x,p,q){
	gamma(p+q)/(gamma(p)*gamma(q))*x^(p-1)*(1-x)^(q-1)
}
x <- seq(0, 1, length.out = n)
mu <- (6 * beta(x, 30, 17) + 4 * beta(x, 3, 11))/10;
mu <- (mu - min(mu))/(max(mu) - min(mu))
noise <- rnorm(n)
y <- mu + sigma * noise

# correlation matrix assumed known and equal to the identity
fit.d <- eBsc(y, method = "D", R0 = diag(n))
plot(fit.d, full = FALSE)
}

\references{
Serra, P. and Krivobokova, T. (2015)\cr
Adaptive Empirical Bayesian Smoothing Splines
}

\keyword{package}
\seealso{
\code{\link[stats]{stl}} (package stats),
\code{\link[stats]{HoltWinters}} (package stats)
}
