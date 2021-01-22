\name{eBsc}
\alias{eBsc}
\title{
Empirical Bayes Smoothing Splines with Correlated Errors
}
\description{
Empirical Bayes smoothing splines with correlated errors. The method uses a recursive algorithm for signal extraction with a non-parametric estimation of the correlation matrix of the errors.
}
\usage{
eBsc(y, q, method, parallel, R0, zero_range, ARpMAq, trace, tol.lambda, tol.rho, max.iter)
}

\arguments{
  \item{y}{Is a univariate numeric vector without missing values. 
}
  \item{q}{Is the value of q if known. If left empty the method considers all possibles q's between 1 and 6 and selects the best one according to the Tq criteria. q=NULL is the default.
}
  \item{method}{Is a method used for the fit. It can take the values "D" (deterministic fit), "P" (parametric fit) and "N" (non-parametric fit).  For example: i) to fit a model with known correlation matrix R.known one should select method = "D" and R0 = R.known; ii) to fit a model with a nonparametric estimation of the correlation and a starting correlation matrix R.start, one should select  method = "N" and R0 = R.start; and iii) to fit a model with an ARMA parametric structure R.ARMA, one should select method="P" and ARpMAq=c(1,0). method = "N" is the default.
}
  \item{parallel}{Is a logical parameter indicating if parallel computation should be used. parallel=FALSE is the default.  	  
}
  \item{R0}{Is the starting correlation matrix. If method = "D" this matrix is not changed by the algorithm.  	  
}
  \item{zero_range}{Is the interval to look for zeros in the estimating equation for the smoothing parameter (lambda).  	  
}
  \item{ARpMAq}{Is the desired ARMA structure for the noise process.  	  
}
  \item{trace}{If true, the process of the algorithm is traced and reported.  	  
}
  \item{tol.lambda}{Tolerance level for lambda.  	  
}
  \item{tol.rho}{Tolerance level for rho.  	  
}
  \item{max.iter}{Maximum number of iterations.  	  
}
}

\details{
The method assumes the data is equidistant. 
}
\value{
A list object of class \code{eBsc} containing the following information.
\item{q.hat}{estimadted q}
\item{lambda.hat}{estimated lambda}
\item{R.hat}{estimated correlation matrix}
\item{f.hat}{estimated function}
\item{f.hat}{estimated variance}
\item{etq.hat}{estimating equation for q}
\item{data}{data used to fit the model}
\item{call}{Call of eBsc}
}
\references{
Serra, P. and Krivobokova, T. (2015)\cr
Adaptive Empirical Bayesian Smoothing Splines
}
\author{
Francisco Rosales, Paulo Serra, Tatyana Krivobokova
}
\seealso{
\code{\link[stats]{stl}} (package stats),
\code{\link[stats]{HoltWinters}} (package stats)
}
\examples{

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
fit <- eBsc(y,  method = "D", R0 = diag(n))
plot(fit, full=FALSE)

}
\keyword{bayesian}
\keyword{non-parametric}
\keyword{correlated errors}