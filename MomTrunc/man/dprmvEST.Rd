\name{dprmvEST}
\alias{dmvEST}
\alias{pmvEST}
\alias{rmvEST}
\title{Multivariate Extended-Skew t Density, Probablilities and Random Deviates Generator}
\description{
  These functions provide the density function, probabilities and a random number
  generator for the multivariate extended-skew t (EST) distribution with mean vector \code{mu}, scale matrix \code{Sigma}, skewness parameter \code{lambda}, extension parameter \code{tau} and degrees of freedom \code{nu}.
}
\usage{
dmvEST(x,mu=rep(0,length(lambda)),Sigma=diag(length(lambda)),lambda,tau=0,nu)
pmvEST(lower = rep(-Inf,length(lambda)),upper=rep(Inf,length(lambda)),
        mu = rep(0,length(lambda)),Sigma,lambda,tau,nu,log2 = FALSE)
rmvEST(n,mu=rep(0,length(lambda)),Sigma=diag(length(lambda)),lambda,tau,nu)
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
  a numeric vector of length \eqn{p} representing the skewness parameter for ST and EST cases. If \code{lambda == 0}, the EST/ST reduces to a t (symmetric) distribution.
  }
  \item{tau}{
  It represents the extension parameter for the EST distribution. If \code{tau == 0}, the EST reduces to a ST distribution.
  }
    \item{nu}{
It represents the degrees of freedom of the Student's t-distribution.
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
\code{\link{dmvST}}, \code{\link{pmvST}}, \code{\link{rmvST}}, \code{\link{meanvarFMD}},\code{\link{meanvarTMD}},\code{\link{momentsTMD}}
}
\examples{
#Univariate case
dmvEST(x = -1,mu = 2,Sigma = 5,lambda = -2,tau = 0.5,nu=4)
rmvEST(n = 100,mu = 2,Sigma = 5,lambda = -2,tau = 0.5,nu=4)
#Multivariate case
mu = c(0.1,0.2,0.3,0.4)
Sigma = matrix(data = c(1,0.2,0.3,0.1,0.2,1,0.4,-0.1,0.3,0.4,1,0.2,0.1,-0.1,0.2,1),
               nrow = length(mu),ncol = length(mu),byrow = TRUE)
lambda = c(-2,0,1,2)
tau = 2
#One observation
dmvEST(x = c(-2,-1,0,1),mu,Sigma,lambda,tau,nu=4)
rmvEST(n = 100,mu,Sigma,lambda,tau,nu=4)
#Many observations as matrix
x = matrix(rnorm(4*10),ncol = 4,byrow = TRUE)
dmvEST(x = x,mu,Sigma,lambda,tau,nu=4)

lower = rep(-Inf,4)
upper = c(-1,0,2,5)
pmvEST(lower,upper,mu,Sigma,lambda,tau,nu=4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Probability}% use one of RShowDoc("KEYWORDS")
\keyword{Multivariate}% __ONLY ONE__ keyword per line
\keyword{t}
\keyword{Student}
\keyword{Skew}
\keyword{Extended}
