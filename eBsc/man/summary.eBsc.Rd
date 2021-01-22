\name{summary.eBsc}
\alias{summary.eBsc}
\title{eBsc Summary}
\description{Takes an \code{eBsc} object produced by \code{eBsc} and summarizes the information of the errors. 
}
\usage{
\method{summary}{eBsc}(object,...)
}
\arguments{
  \item{object}{\code{eBsc} object. 
  }
  \item{\dots}{further arguments to be passed to summary().
} 
}
\value{
The function gives basic statistics of the error from applying \code{eBsc}. 
}
\references{
Serra, P. and Krivobokova, T. (2015)\cr
Adaptive Empirical Bayesian Smoothing Splines
}
\author{
Francisco Rosales, Paulo Serra, Tatyana Krivobokova
}

\seealso{
\code{\link[eBsc]{plot.eBsc}} (package eBsc),
}
\examples{
# simulated data	
library(eBsc)
n <- 250
sigma <- 0.05

Basis <- list()
for(i in 1:6){Basis[[i]] <- drbasis(nn = n, qq = i)}
coef3 <- c(rep(0,3),(pi*(2:(n-2)))^(-3.1)) * (cos(2*(1:n)))
A3 <- Basis[[3]]$eigenvectors
mu <- - A3\%*\%coef3 
mu <- (mu - min(mu))/(max(mu) - min(mu))
noise <- rnorm(n)
y <- mu + sigma * noise

# correlation matrix assumed known and equal to the identity
fit <- eBsc(y, method = "D", R0 = diag(n))

summary(fit)
}
\keyword{summary}
\keyword{correlated errors}
