\name{pmvnormt}
\alias{pmvnormt}
\alias{pmvnormt}
\title{Multivariate normal and Student-t probabilities}
\description{
  Computation of Multivariate normal and Student-t probabilities using the classic Genz method form packages \code{mvtnorm} and \code{tlrmvnmvt} packages. In order to save computational effort, it chooses whether to use the function \code{pmvtnorm} (\code{pmvt}) from \code{mvtnorm}, or functions \code{pmvn} (\code{pmvt}) from the \code{tlrmvnmvt} package, depending of the vector size \code{p}, real or integer degrees of freedom \code{nu}.
}
\usage{
pmvnormt(lower = rep(-Inf,ncol(sigma)),upper = rep(Inf,ncol(sigma)),
mean = rep(0,ncol(sigma)),sigma,nu = NULL,uselog2 = FALSE)
}
\arguments{
  \item{lower}{lower integration limits, a numeric vector of length p}
  \item{upper}{upper integration limits, a numeric vector of length p}
  \item{mean}{the mean parameter, a numeric vector of length p}
  \item{sigma}{the covariance matrix, a square matrix that matches the length of `lower`}
  \item{nu}{degrees of freedom, a positive real number. If \code{NULL}, normal case is considered}
  \item{uselog2}{a boolean variable, indicating if the log2 result should be returned.
    This is useful when the true probability is too small for the machine precision}
}
\value{
The estimated probability or its log2 if \code{uselog2 == TRUE}
}


\references{
  Genz, A. (1992), "Numerical computation of multivariate
normal probabilities," Journal of Computational and
Graphical Statistics, 1, 141-149.

Cao, J., Genton, M. G., Keyes, D. E., & Turkiyyah, G. M. "Exploiting Low Rank Covariance Structures for Computing High-Dimensional Normal and Student- t Probabilities" (2019) <https://marcgenton.github.io/2019.CGKT.manuscript.pdf>
}


\author{Christian E. Galarza <\email{cgalarza88@gmail.com}> and
Victor H. Lachos <\email{hlachos@uconn.edu}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}
\note{
If \code{is.null(nu)}, normal case is considered.
}

\seealso{
\code{\link{onlymeanTMD}},\code{\link{meanvarTMD}},\code{\link{momentsFMD}},\code{\link{momentsTMD}},\code{\link{meanvarFMD}},\code{\link{dmvSN}},\code{\link{pmvSN}},\code{\link{rmvSN}}, \code{\link{dmvESN}},\code{\link{pmvESN}},\code{\link{rmvESN}}, \code{\link{dmvST}},\code{\link{pmvST}},\code{\link{rmvST}}, \code{\link{dmvEST}},\code{\link{pmvEST}},\code{\link{rmvEST}}
}

\examples{
a = c(-0.8,-0.7,-0.6)
b = c(0.5,0.6,0.7)
mu = c(0.1,0.2,0.3)
Sigma = matrix(data = c(1,0.2,0.3,0.2,1,0.4,0.3,0.4,1),
               nrow = length(mu),ncol = length(mu),byrow = TRUE)

pmvnormt(lower = a,upper = b,mean = mu,sigma = Sigma) #normal case
pmvnormt(lower = a,upper = b,mean = mu,sigma = Sigma,nu = 4.23)  #t case 
pmvnormt(lower = a,upper = b,mean = mu,sigma = Sigma,nu = 4.23,uselog2 = TRUE) 

}