\name{dprmvESN}
\alias{dmvESN}
\alias{pmvESN}
\alias{rmvESN}
\title{Multivariate Extended-Skew Normal Density, Probablilities and Random Deviates Generator}
\description{
  These functions provide the density function, probabilities and a random number
  generator for the multivariate extended-skew normal (ESN) distribution with mean vector \code{mu}, scale matrix \code{Sigma}, skewness parameter \code{lambda} and extension parameter \code{tau}.
}
\usage{
dmvESN(x,mu=rep(0,length(lambda)),Sigma=diag(length(lambda)),lambda,tau=0)
pmvESN(lower = rep(-Inf,length(lambda)),upper=rep(Inf,length(lambda)),
        mu = rep(0,length(lambda)),Sigma,lambda,tau,log2 = FALSE)
rmvESN(n,mu=rep(0,length(lambda)),Sigma=diag(length(lambda)),lambda,tau=0)
}
\arguments{
   \item{x}{vector or matrix of quantiles. If \code{x} is a matrix, each
     row is taken to be a quantile.}
   \item{n}{number of observations.}
   \item{lower}{
the vector of lower limits of length \eqn{p}.
}
\item{upper}{
the vector of upper limits of length \eqn{p}.
}
   \item{mu}{
  a numeric vector of length \eqn{p} representing the location parameter.
  }
  \item{Sigma}{
  a numeric positive definite matrix with dimension \eqn{p}x\eqn{p} representing the scale parameter.
  }
  \item{lambda}{
  a numeric vector of length \eqn{p} representing the skewness parameter for SN and ESN cases. If \code{lambda == 0}, the ESN/SN reduces to a normal (symmetric) distribution.
  }
  \item{tau}{
  It represents the extension parameter for the ESN distribution. If \code{tau == 0}, the ESN reduces to a SN distribution.
  }
    \item{log2}{
  a boolean variable, indicating if the log2 result should be returned. This is useful when the true probability is too small for the machine precision.
  }
}

\references{
Cao, J., Genton, M. G., Keyes, D. E., & Turkiyyah, G. M. (2019) "Exploiting Low
	Rank Covariance Structures for Computing High-Dimensional Normal and Student-
	t Probabilities" <\url{https://marcgenton.github.io/2019.CGKT.manuscript.pdf}>.
	
Galarza C.E., Matos L.A., Dey D.K. & Lachos V.H. (2019) On Moments of Folded and Truncated Multivariate Extended Skew-Normal Distributions. Technical report. ID 19-14. University of Connecticut <\url{https://stat.uconn.edu/tech-reports-2019/}>.

Genz, A., "Numerical computation of multivariate normal probabilities," Journal of Computational and Graphical Statistics, 1, 141-149 (1992) <doi:10.1080/10618600.1992.10477010>.
}

\author{Christian E. Galarza <\email{cgalarza88@gmail.com}> and
Victor H. Lachos <\email{hlachos@uconn.edu}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{dmvSN}}, \code{\link{pmvSN}}, \code{\link{rmvSN}}, \code{\link{meanvarFMD}},\code{\link{meanvarTMD}},\code{\link{momentsTMD}}
}
\examples{
#Univariate case
dmvESN(x = -1,mu = 2,Sigma = 5,lambda = -2,tau = 0.5)
rmvESN(n = 100,mu = 2,Sigma = 5,lambda = -2,tau = 0.5)
#Multivariate case
mu = c(0.1,0.2,0.3,0.4)
Sigma = matrix(data = c(1,0.2,0.3,0.1,0.2,1,0.4,-0.1,0.3,0.4,1,0.2,0.1,-0.1,0.2,1),
               nrow = length(mu),ncol = length(mu),byrow = TRUE)
lambda = c(-2,0,1,2)
tau = 2
#One observation
dmvESN(x = c(-2,-1,0,1),mu,Sigma,lambda,tau)
rmvESN(n = 100,mu,Sigma,lambda,tau)
#Many observations as matrix
x = matrix(rnorm(4*10),ncol = 4,byrow = TRUE)
dmvESN(x = x,mu,Sigma,lambda,tau)

lower = rep(-Inf,4)
upper = c(-1,0,2,5)
pmvESN(lower,upper,mu,Sigma,lambda,tau)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Probability}% use one of RShowDoc("KEYWORDS")
\keyword{Multivariate}% __ONLY ONE__ keyword per line
\keyword{Normal}
\keyword{Skew}
\keyword{Extended}
