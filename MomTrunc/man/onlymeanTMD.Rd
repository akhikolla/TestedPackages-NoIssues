\name{onlymeanTMD}
\alias{onlymeanTMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean for doubly truncated multivariate distributions
}
\description{
It computes the mean vector for some doubly truncated skew-elliptical distributions. It supports the \code{p}-variate Normal, Skew-normal (SN), Extended Skew-normal (ESN) and Unified Skew-normal (SUN) as well as the Student's-t, Skew-t (ST), Extended Skew-t (EST) and Unified Skew-t (SUT) distribution.
}
\usage{
onlymeanTMD(lower = rep(-Inf, length(mu)),upper = rep(Inf,length(mu)),mu,Sigma,
             lambda = NULL,tau = NULL,Gamma = NULL,nu = NULL,dist)
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
a numeric vector of length \eqn{p} representing the skewness parameter for SN and ESN cases. If \code{lambda == 0}, the ESN/SN reduces to a normal (symmetric) distribution.
}
  \item{tau}{
It represents the extension parameter for the ESN distribution. If \code{tau == 0}, the ESN reduces to a SN distribution.
}
 \item{Gamma}{
a correlation matrix with dimension \eqn{q}x\eqn{q}. It must be provided only for the SUN and SUT cases. For particular cases SN, ESN, ST and EST, we have that \code{Gamma == 1} (see examples at the end of this help).
}
  \item{nu}{
It represents the degrees of freedom for the Student's t-distribution.
}
  \item{dist}{
represents the truncated distribution to be used. The values are \code{normal}, \code{SN} , \code{ESN} and \code{SUN} for the doubly truncated Normal, Skew-normal, Extended Skew-normal and Unified-skew normal distributions and,
\code{t}, \code{ST} , \code{EST} and \code{SUT} for the for the doubly truncated Student-t, Skew-t, Extended Skew-t and Unified skew-t distributions.
}
}
\details{
Univariate case is also considered, where \code{Sigma} will be the variance \eqn{\sigma^2}. Normal case code is an R adaptation of the Matlab available function \code{dtmvnmom.m} from Kan & Robotti (2017) and it is used for \code{p<=3}. For higher dimensions we use an extension of the algorithm in Vaida (2009).
}
\value{
It returns the mean vector of length \eqn{p}.
}
\references{
Arellano-Valle, R. B. & Genton, M. G. (2005). On fundamental skew distributions. Journal of Multivariate
Analysis, 96, 93-116.

Arellano-Valle, R. B. & Genton, M. G. (2010). Multivariate extended skew-t distributions and related
families. Metron, 68(3), 201-234.

Galarza C.E., Matos L.A., Dey D.K. & Lachos V.H. (2019) On Moments of Folded and Truncated Multivariate Extended Skew-Normal Distributions. Technical report. ID 19-14. University of Connecticut.

Kan, R., & Robotti, C. (2017). On moments of folded and truncated multivariate normal distributions. Journal of Computational and Graphical Statistics, 26(4), 930-934.

Kirkby J. Lars, Nguyen D. and Nguyen D. (2019). Moments of Student's t-distribution: A Unified Approach. <\url{https://arxiv.org/abs/1912.01607}>

Vaida, F. & Liu, L. (2009). Fast implementation for normal mixed effects models with censored response. Journal of Computational and Graphical Statistics, 18(4), 797-817.
}
\author{Christian E. Galarza <\email{cgalarza88@gmail.com}> and
Victor H. Lachos <\email{hlachos@uconn.edu}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}
\note{
Degrees of freedom must be a positive integer. If \code{nu >= 300}, Normal case is considered."
}

\seealso{
\code{\link{momentsTMD}}, \code{\link{meanvarFMD}}, \code{\link{momentsFMD}},\code{\link{dmvESN}},\code{\link{rmvESN}}
}
\examples{
a = c(-0.8,-0.7,-0.6)
b = c(0.5,0.6,0.7)
mu = c(0.1,0.2,0.3)
Sigma = matrix(data = c(1,0.2,0.3,0.2,1,0.4,0.3,0.4,1),
               nrow = length(mu),ncol = length(mu),byrow = TRUE)
value1 = onlymeanTMD(a,b,mu,Sigma,dist="normal")

# Now works for for any nu>0
value2 = onlymeanTMD(a,b,mu,Sigma,dist = "t",nu = 0.87)

value3 = onlymeanTMD(a,b,mu,Sigma,lambda = c(-2,0,1),dist = "SN")
value4 = onlymeanTMD(a,b,mu,Sigma,lambda = c(-2,0,1),tau = 1,dist = "ESN")
value5 = onlymeanTMD(a,b,mu,Sigma,lambda = c(-2,0,1),tau = 1,nu = 4,dist = "EST")

#Skew-unified Normal (SUN) and Skew-unified t (SUT) distributions

Lambda = matrix(c(1,0,2,-3,0,-1),3,2) #A skewness matrix p times q
Gamma  = matrix(c(1,-0.5,-0.5,1),2,2) #A correlation matrix q times q
tau    = c(-1,2)                      #A vector of extension parameters of dim q

value6 = onlymeanTMD(a,b,mu,Sigma,lambda = Lambda,tau = c(-1,2),Gamma = Gamma,dist = "SUN")
value7 = onlymeanTMD(a,b,mu,Sigma,lambda = Lambda,tau = c(-1,2),Gamma = Gamma,nu = 4,dist = "SUT")

#The ESN and EST as particular cases of the SUN and SUT for q=1

Lambda = matrix(c(-2,0,1),3,1)
Gamma  = 1
value8 = onlymeanTMD(a,b,mu,Sigma,lambda = Lambda,tau = 1,Gamma = Gamma,dist = "SUN")
value9 = onlymeanTMD(a,b,mu,Sigma,lambda = Lambda,tau = 1,Gamma = Gamma,nu = 4,dist = "SUT")

round(value4,2) == round(value8,2)
round(value5,2) == round(value9,2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Truncated}% use one of  RShowDoc("KEYWORDS")
\keyword{Multivariate}% __ONLY ONE__ keyword per line
\keyword{Student's t}
\keyword{Normal}
\keyword{Skew}
\keyword{Extended}
