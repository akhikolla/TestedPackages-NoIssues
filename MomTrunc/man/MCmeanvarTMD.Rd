\name{MCmeanvarTMD}
\alias{MCmeanvarTMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Monte Carlo Mean and variance for doubly truncated multivariate distributions
}
\description{
  It computes the Monte Carlo mean vector and variance-covariance matrix for some doubly truncated skew-elliptical distributions. Monte Carlo simulations are performed via slice Sampling.
  It supports the \code{p}-variate Normal, Skew-normal (SN), Extended Skew-normal (ESN) and Unified Skew-normal (SUN) as well as the Student's-t, Skew-t (ST), Extended Skew-t (EST) and Unified Skew-t (SUT) distribution.
}
\usage{
MCmeanvarTMD(lower = rep(-Inf,length(mu)),upper = rep(Inf,length(mu)),mu,Sigma
,lambda = NULL,tau = NULL,Gamma = NULL,nu = NULL,dist,n = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
a numeric matrix of dimension \eqn{p}x\eqn{q} representing the skewness/shape matrix parameter for the SUN and SUT distribution. For the ESN and EST distributions (\eqn{q=1}), \code{lambda} is a numeric vector of dimension \eqn{p} (see examples at the end of this help). If \code{all(lambda == 0)}, the SUN/ESN/SN (SUT/EST/ST) reduces to a normal (t) symmetric distribution.
}
  \item{tau}{
a numeric vector of length \eqn{q} representing the extension parameter for the SUN and SUT distribution. For the ESN and EST distributions, \code{tau} is a positive scalar (\eqn{q=1}). Furthermore, if \code{tau == 0}, the ESN (EST) reduces to a SN (ST) distribution.
}
  \item{Gamma}{
a correlation matrix with dimension \eqn{q}x\eqn{q}. It must be provided only for the SUN and SUT cases. For particular cases SN, ESN, ST and EST, we have that \code{Gamma == 1} (see examples at the end of this help).
}
  \item{nu}{
It represents the degrees of freedom for the Student's t-distribution being a positive real number.
}
\item{dist}{
  represents the truncated distribution to be used. The values are \code{normal}, \code{SN} , \code{ESN} and \code{SUN} for the doubly truncated Normal, Skew-normal, Extended Skew-normal and Unified-skew normal distributions and,
  \code{t}, \code{ST} , \code{EST} and \code{SUT} for the for the doubly truncated Student-t, Skew-t, Extended Skew-t and Unified skew-t distributions.
}
\item{n}{
  number of Monte Carlo samples to be generated.
}

}
\value{
  It returns a list with three elements:
    \item{mean}{the estimate for the mean vector of length \eqn{p}}
  \item{EYY}{the estimate for the second moment matrix of dimensions \eqn{p}x\eqn{p}}
  \item{varcov}{the estimate for the variance-covariance matrix of dimensions \eqn{p}x\eqn{p}}
}
\references{
  Arellano-Valle, R. B. & Genton, M. G. (2005). On fundamental skew distributions. Journal of Multivariate
  Analysis, 96, 93-116.
  
  Ho, H. J., Lin, T. I., Chen, H. Y., & Wang, W. L. (2012). Some results on the truncated multivariate t distribution. Journal of Statistical Planning and Inference, 142(1), 25-40.
}

\author{Christian E. Galarza <\email{cgalarza88@gmail.com}> and
Victor H. Lachos <\email{hlachos@uconn.edu}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}

\seealso{
\code{\link{meanvarTMD}}, \code{\link{rmvSN}},\code{\link{rmvESN}},\code{\link{rmvST}}, \code{\link{rmvEST}}
}
\examples{
a = c(-0.8,-0.7,-0.6)
b = c(0.5,0.6,0.7)
mu = c(0.1,0.2,0.3)
Sigma = matrix(data = c(1,0.2,0.3,0.2,1,0.4,0.3,0.4,1),
               nrow = length(mu),ncol = length(mu),byrow = TRUE)


## Normal case

# Theoretical value
value1 = meanvarTMD(a,b,mu,Sigma,dist="normal")

#MC estimate
MC11 = MCmeanvarTMD(a,b,mu,Sigma,dist="normal") #by defalut n = 10000
MC12 = MCmeanvarTMD(a,b,mu,Sigma,dist="normal",n = 10^5) #more precision

## Skew-t case
 
# Theoretical value
value2 = meanvarTMD(a,b,mu,Sigma,lambda = c(-2,0,1),nu = 4,dist = "ST")

#MC estimate
MC21 = MCmeanvarTMD(a,b,mu,Sigma,lambda = c(-2,0,1),nu = 4,dist = "ST")

## More...

MC5 = MCmeanvarTMD(a,b,mu,Sigma,lambda = c(-2,0,1),tau = 1,dist = "ESN")
MC6 = MCmeanvarTMD(a,b,mu,Sigma,lambda = c(-2,0,1),tau = 1,nu = 4,dist = "EST")

\dontrun{
#Skew-unified Normal (SUN) and Skew-unified t (SUT) distributions

Lambda = matrix(c(1,0,2,-3,0,-1),3,2) #A skewness matrix p times q
Gamma  = matrix(c(1,-0.5,-0.5,1),2,2) #A correlation matrix q times q
tau    = c(-1,2)                      #A vector of extension parameters of dim q

MC7 = MCmeanvarTMD(a,b,mu,Sigma,lambda = Lambda,tau = c(-1,2),Gamma = Gamma,dist = "SUN")
MC8 = MCmeanvarTMD(a,b,mu,Sigma,lambda = Lambda,tau = c(-1,2),Gamma = Gamma,nu = 1,dist = "SUT")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Monte Carlo}
\keyword{Truncated}% use one of  RShowDoc("KEYWORDS")
\keyword{Multivariate}% __ONLY ONE__ keyword per line
\keyword{Student's t}
\keyword{Normal}
\keyword{Skew}
\keyword{Extended}
\keyword{Unified}
\keyword{Selection}