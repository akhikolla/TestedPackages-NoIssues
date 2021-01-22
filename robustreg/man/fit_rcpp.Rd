\name{fit_rcpp}
\alias{fit_rcpp}
\docType{package}
\title{
Predict y from X and b
}
\description{
Predict y vector from X design matrix and b vector 
}
\usage{
fit_rcpp(X,b)
}
\arguments{
\item{X}{Design matrix}
\item{b}{Estimates of beta}
}
\author{Ian M. Johnson}
\examples{
j <- rep(1, 5)
x1 <- rnorm(5)
x2 <- rnorm(5, 10, 20)
X = as.matrix(data.frame(j, x1, x2))
b <- 1:3
fit_rcpp(X, b)
}
