\name{plot.eBsc}
\alias{plot.eBsc}
\title{Plot fitted components
}
\description{
Plot fitted components and the acf of the errors.
}
\usage{
\method{plot}{eBsc}(x,full=FALSE,...)
}
\arguments{
  \item{x}{\code{eBsc}  object.
}
  \item{full}{plot option. If TRUE graphial details of the estimation are provided. If FALSE a simple plot of the estimation and its confidence bands is provided. 
}
  \item{\dots}{further arguments to be passed to plot().
}

}
\details{
if the \code{eBsc} plots the fits and the acf of the errors.
}
\value{
The function returns the selected plots. 
}
\references{
Serra, P. and Krivobokova, T. (2015)\cr
Adaptive Empirical Bayesian Smoothing Splines
}
\author{
Francisco Rosales, Paulo Serra, Tatyana Krivobokova. 
}

\examples{

library(eBsc)
n <- 250
sigma <- 0.05
Basis <- list()
for(i in 1:6) Basis[[i]] <- drbasis(nn = n, qq = i)
coef3 <- c(rep(0,3),(pi*(2:(n-2)))^(-3.1))*(cos(2*(1:n)))
A3 <- Basis[[3]]$eigenvectors
mu <- A3\%*\%coef3
mu <- (mu-min(mu))/(max(mu)-min(mu))
noise <- rnorm(n)
y <- mu + sigma * noise

#correlation assumed known and equal to the identity
fit.d <- eBsc(y, method = "D", R0 = diag(n))

#simple plot by  
plot(fit.d, full = FALSE)

}

\keyword{plot}

