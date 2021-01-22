\name{meanvarTMD}
\alias{meanvarTMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean and variance for doubly truncated multivariate distributions
}
\description{
It computes the mean vector and variance-covariance matrix for some doubly truncated skew-elliptical distributions. It supports the \code{p}-variate Normal, Skew-normal (SN), Extended Skew-normal (ESN) and Unified Skew-normal (SUN) as well as the Student's-t, Skew-t (ST), Extended Skew-t (EST) and Unified Skew-t (SUT) distribution.
}
\usage{
meanvarTMD(lower = rep(-Inf,length(mu)),upper = rep(Inf,length(mu)),mu,Sigma
,lambda = NULL,tau = NULL,Gamma = NULL,nu = NULL,dist)
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
}
\details{
Univariate case is also considered, where \code{Sigma} will be the variance \eqn{\sigma^2}. Normal case code is an R adaptation of the Matlab available function \code{dtmvnmom.m} from Kan & Robotti (2017) and it is used for \code{p<=3}. For higher dimensions we use an extension of the algorithm in Vaida (2009).
}
\value{
It returns a list with three elements:
\item{mean}{the mean vector of length \eqn{p}}
\item{EYY}{the second moment matrix of dimensions \eqn{p}x\eqn{p}}
\item{varcov}{the variance-covariance matrix of dimensions \eqn{p}x\eqn{p}}
}
\references{
Arellano-Valle, R. B. & Genton, M. G. (2005). On fundamental skew distributions. Journal of Multivariate
Analysis, 96, 93-116.

Arellano-Valle, R. B. & Genton, M. G. (2010). Multivariate extended skew-t distributions and related
families. Metron, 68(3), 201-234.

Galarza C.E., Matos L.A., Dey D.K. & Lachos V.H. (2019) On Moments of Folded and Truncated Multivariate Extended Skew-Normal Distributions. Technical report. ID 19-14. University of Connecticut <\url{https://stat.uconn.edu/tech-reports-2019/}>.

Ho, H. J., Lin, T. I., Chen, H. Y., & Wang, W. L. (2012). Some results on the truncated multivariate t distribution. Journal of Statistical Planning and Inference, 142(1), 25-40.

Kan, R., & Robotti, C. (2017). On moments of folded and truncated multivariate normal distributions. Journal of Computational and Graphical Statistics, 26(4), 930-934.

Kirkby J. Lars, Nguyen D. and Nguyen D. (2019). Moments of Student's t-distribution: A Unified Approach. <\url{https://arxiv.org/abs/1912.01607}>

Vaida, F. & Liu, L. (2009). Fast implementation for normal mixed effects models with censored response. Journal of Computational and Graphical Statistics, 18(4), 797-817.
}
\author{Christian E. Galarza <\email{cgalarza88@gmail.com}> and
Victor H. Lachos <\email{hlachos@uconn.edu}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}
\note{
If \code{nu >= 300}, Normal case is considered."
}

\section{Warning }{For the \eqn{t} cases, the algorithm supports degrees of freedom \code{nu <= 2}.}

\seealso{
\code{\link{MCmeanvarTMD}},
\code{\link{momentsTMD}}, \code{\link{meanvarFMD}}, \code{\link{meanvarFMD}},\code{\link{momentsFMD}}, \code{\link{dmvSN}},\code{\link{pmvSN}},\code{\link{rmvSN}}, \code{\link{dmvESN}},\code{\link{pmvESN}},\code{\link{rmvESN}}, \code{\link{dmvST}},\code{\link{pmvST}},\code{\link{rmvST}}, \code{\link{dmvEST}},\code{\link{pmvEST}},\code{\link{rmvEST}}
}
\examples{
a = c(-0.8,-0.7,-0.6)
b = c(0.5,0.6,0.7)
mu = c(0.1,0.2,0.3)
Sigma = matrix(data = c(1,0.2,0.3,0.2,1,0.4,0.3,0.4,1),
               nrow = length(mu),ncol = length(mu),byrow = TRUE)

# Theoretical value
value1 = meanvarTMD(a,b,mu,Sigma,dist="normal")

#MC estimate
MC11 = MCmeanvarTMD(a,b,mu,Sigma,dist="normal") #by defalut n = 10000
MC12 = MCmeanvarTMD(a,b,mu,Sigma,dist="normal",n = 10^5) #more precision

# Now works for for any nu>0
value2 = meanvarTMD(a,b,mu,Sigma,dist = "t",nu = 0.87)

value3 = meanvarTMD(a,b,mu,Sigma,lambda = c(-2,0,1),dist = "SN")
value4 = meanvarTMD(a,b,mu,Sigma,lambda = c(-2,0,1),nu = 4,dist = "ST")
value5 = meanvarTMD(a,b,mu,Sigma,lambda = c(-2,0,1),tau = 1,dist = "ESN")
value6 = meanvarTMD(a,b,mu,Sigma,lambda = c(-2,0,1),tau = 1,nu = 4,dist = "EST")

\dontrun{
#Skew-unified Normal (SUN) and Skew-unified t (SUT) distributions

Lambda = matrix(c(1,0,2,-3,0,-1),3,2) #A skewness matrix p times q
Gamma  = matrix(c(1,-0.5,-0.5,1),2,2) #A correlation matrix q times q
tau    = c(-1,2)                      #A vector of extension parameters of dim q

value7 = meanvarTMD(a,b,mu,Sigma,lambda = Lambda,tau = c(-1,2),Gamma = Gamma,dist = "SUN")
value8 = meanvarTMD(a,b,mu,Sigma,lambda = Lambda,tau = c(-1,2),Gamma = Gamma,nu = 4,dist = "SUT")


#The ESN and EST as particular cases of the SUN and SUT for q=1

Lambda = matrix(c(-2,0,1),3,1)
Gamma  = 1
value9 = meanvarTMD(a,b,mu,Sigma,lambda = Lambda,tau = 1,Gamma = Gamma,dist = "SUN")
value10 = meanvarTMD(a,b,mu,Sigma,lambda = Lambda,tau = 1,Gamma = Gamma,nu = 4,dist = "SUT")

round(value5$varcov,2) == round(value9$varcov,2)
round(value6$varcov,2) == round(value10$varcov,2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Truncated}% use one of  RShowDoc("KEYWORDS")
\keyword{Multivariate}% __ONLY ONE__ keyword per line
\keyword{Student's t}
\keyword{Normal}
\keyword{Skew}
\keyword{Extended}
\keyword{Unified}
\keyword{Selection}