\name{drbasis}
\alias{drbasis}
\title{
Computation of the Demmler-Reinsch basis.
}
\description{
The Demmler-Reinsch basis is provided for a given smoothness degree in a uniform grid. 
}
\usage{
drbasis(nn,qq)
}
\arguments{
  \item{nn}{Number of design points in the uniform grid. 
}
  \item{qq}{Smoothness degree of the basis.
}
}
\details{
The use of large numbers required by the basis is handled by the package Brobdingnag. The method assumes the grid is equidistant. Missing values are not supported.
}
\value{
A list object containing the following information.
\item{eigenvalues}{estimadted eigenvalues}
\item{eigenvectors}{estimated eigenvectors}
\item{eigenvectorsQR}{orthonormal eigenvectors}
\item{x}{equidistant grid used to build the basis}
}
\references{
Rosales F. (2016). \cr
Empirical Bayesian Smoothing Splines for Signals with Correlated Errors: Methods and Applications

Serra, P. and Krivobokova, T. (2015)\cr
Adaptive Empirical Bayesian Smoothing Splines
}
\author{
Francisco Rosales
}
\examples{
oldpar <- par(no.readonly = TRUE)
#plot elements of the basis
library(eBsc)
n <- 100
Basis <- list()
for(i in 1:6){Basis[[i]] <- drbasis( nn = n, qq = i)}

#eigenvalues
par(mfrow = c(3,2), mar = c(4,2,2,2))
for(i in 1:6){
name <- paste("Eigenvalues (q = ",i,")", sep = "")
plot(Basis[[i]]$eigenvalues,
type = 'l', lwd = 2, xlab = "x", ylab = "", main = name)	
}
par(oldpar) 

#eigenvectors for q = 3
par(mfrow = c(3,2), mar = c(4,2,2,2))
for(i in 1:6){
name <- paste("Eigenvector n. ", i + 3, sep = "")
plot(Basis[[i]]$eigenvectorsQR[, i + 3],
type = 'l', lwd = 2, xlab = "x", ylab = "", main = name)	
}
par(oldpar) 

#example of a smooth function in the Demmler-Reinsch basis
library(eBsc)
n <- 200
Basis <- list()
for(i in 1:6){Basis[[i]] <- drbasis(nn = n, qq = i)}
coef3 <- c(rep(0,3), (pi*(2:(n-2))) ^ (-3.1)) * (cos(2*(1:n)))
A3 <- Basis[[3]]$eigenvectors
mu <- -A3\%*\%coef3
mu <- (mu - min(mu)) / (max(mu) - min(mu))
plot(mu, xlab = "x", ylab = "mu", type = 'l', lwd = 2)
par(oldpar)

}

\keyword{Demmler-Reinsch basis}
\keyword{non-parametric}