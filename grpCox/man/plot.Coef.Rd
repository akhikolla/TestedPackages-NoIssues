\name{plot.Coef}
\alias{plot.Coef}
\title{
Plots the coefficient paths
}
\description{
Plots the coefficient values as a function of the lambda values used.
}
\usage{
\method{plot}{Coef}(x, lambda, label=TRUE,xlab="log(Lambda)", 
ylab="Coefficients",title=NULL,...)
}
\arguments{
  \item{x}{A matrix of coefficients.}
  \item{lambda}{The lambda values used.}
  \item{label}{The indices of covariates. Default is TRUE.}
  \item{xlab}{The name of the x-axis.}
  \item{ylab}{The name of the y-axis.}
  \item{title}{The title of the plot.}
  \item{\dots}{further arguments to plot}
}
\value{No return value.}
\details{
A plot is produced, and nothing is returned.
}

\author{
Xuan Dang <\email{xuandang11289@gmail.com}>
}

\examples{
set.seed(100001)
N <- 50
p <- 6
times <- 1:p
rho <- 0.5
H <- abs(outer(times, times, "-"))
C <- 1 * rho^H
C[cbind(1:p, 1:p)] <- C[cbind(1:p, 1:p)] 
sigma <- matrix(C,p,p)
mu <- rep(0,p)
x <- mvrnorm(n=N, mu, sigma)

beta <- c(0, .8, 1, 2, 1, 0)
hx <- exp(x \%*\% beta) 
ty <- rexp(N,hx) 
tcens <- 1 - rbinom(n=N, prob = 0.2, size = 1)
y <- data.frame(illt=ty, ills=tcens)
names(y) <- c("time", "status")

group <- list(g1 = c(1,2,3,4), g2 = c(1,2,6), g3 = c(2,3), g4 = c(4,5), g5 = c(5))
fit <- grpCoxOverlap(x, y, group, penalty="glasso", nlambda=50)
# plot the coefficient values in latent space
plot.gCoef(fit$aBetaLatent, fit$glatent, fit$lambda)
# plot the coefficient values in original space
plot.Coef(fit$aBetaOri, fit$lambda)
}

\keyword{Cox models}
\keyword{group regularization}
\keyword{overlapping group}