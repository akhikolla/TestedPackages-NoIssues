\name{plot.llCV}
\alias{plot.llCV}
\title{
Plot the cross-validation curve produced by cv.grpCox or cv.grpCoxOverlap
}
\description{
Plots the cross-validation curve, and upper and lower standard deviation curves, as a function of the lambda values used.
}
\usage{
\method{plot}{llCV}(x,...)
}
\arguments{
  \item{x}{fitted cv.grpCox or cv.grpCoxOverlap object}
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
set.seed(200)
N <- 50
p <- 9
x <- matrix(rnorm(N * p), nrow = N)
beta <- c(.65,.65,0,0,.65,.65,0,.65,0)
hx <- exp(x \%*\% beta) 
ty <- rexp(N,hx) 
tcens <- 1 - rbinom(n=N, prob = 0.2, size = 1)
y <- data.frame(illt=ty, ills=tcens)
names(y) <- c("time", "status")

g <- c(1,1,2,2,3,3,2,3,2)
m <- c(sqrt(2),sqrt(4),sqrt(3))

cvfit <- cv.grpCox(x,y,g,m,penalty="glasso")
plot.llCV(cvfit)
}

\keyword{Cox models}
\keyword{group regularization}
\keyword{overlapping group}
