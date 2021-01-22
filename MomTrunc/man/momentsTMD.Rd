\name{momentsTMD}
\alias{momentsTMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Moments for doubly truncated multivariate distributions
}
\description{
It computes kappa-th order moments for for some doubly truncated skew-elliptical distributions. It supports the \code{p}-variate Normal, Skew-normal (SN) and Extended Skew-normal (ESN), as well as the Student's-t, Skew-t (ST) and the Extended Skew-t (EST) distribution.
}
\usage{
momentsTMD(kappa,lower = rep(-Inf,length(mu)),upper = rep(Inf,length(mu)),mu,Sigma,
lambda = NULL,tau = NULL,nu = NULL,dist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kappa}{
moments vector of length \eqn{p}. All its elements must be integers greater or equal to \eqn{0}. For the Student's-t case, \code{kappa} can be a scalar representing the order of the moment.
}
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
  \item{nu}{
It represents the degrees of freedom for the Student's t-distribution being a positive real number.
}
  \item{dist}{
represents the truncated distribution to be used. The values are \code{normal}, \code{SN} and \code{ESN} for the doubly truncated Normal, Skew-normal and Extended Skew-normal distributions and,
\code{t}, \code{ST} and \code{EST} for the for the doubly truncated Student-t, Skew-t and Extended Skew-t distributions.
}
}
\details{
Univariate case is also considered, where \code{Sigma} will be the variance \eqn{\sigma^2}.
}
\value{
A data frame containing \eqn{p+1} columns. The \eqn{p} first containing the set of combinations of exponents summing up to \code{kappa} and the last column containing the the expected value. Normal cases (ESN, SN and normal) return \code{prod(kappa)+1} moments while the Student's t-distribution case returns all moments of order up to \code{kappa}. See example section.
}

\references{
Arellano-Valle, R. B. & Genton, M. G. (2005). On fundamental skew distributions. Journal of Multivariate
Analysis, 96, 93-116.

Arellano-Valle, R. B. & Genton, M. G. (2010). Multivariate extended skew-t distributions and related
families. Metron, 68(3), 201-234.

Galarza C.E., Matos L.A., Dey D.K. & Lachos V.H. (2019) On Moments of Folded and Truncated Multivariate Extended Skew-Normal Distributions. Technical report. ID 19-14. University of Connecticut <\url{https://stat.uconn.edu/tech-reports-2019/}>.

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

\seealso{
\code{\link{onlymeanTMD}},\code{\link{meanvarTMD}},\code{\link{momentsFMD}},\code{\link{meanvarFMD}},\code{\link{dmvSN}},\code{\link{pmvSN}},\code{\link{rmvSN}}, \code{\link{dmvESN}},\code{\link{pmvESN}},\code{\link{rmvESN}}, \code{\link{dmvST}},\code{\link{pmvST}},\code{\link{rmvST}}, \code{\link{dmvEST}},\code{\link{pmvEST}},\code{\link{rmvEST}}
}

\examples{
a = c(-0.8,-0.7,-0.6)
b = c(0.5,0.6,0.7)
mu = c(0.1,0.2,0.3)
Sigma = matrix(data = c(1,0.2,0.3,0.2,1,0.4,0.3,0.4,1),
               nrow = length(mu),ncol = length(mu),byrow = TRUE)
value1 = momentsTMD(c(2,0,1),a,b,mu,Sigma,dist="normal")
value2 = momentsTMD(c(2,0,1),a,b,mu,Sigma,dist = "t",nu = 7)
value3 = momentsTMD(c(2,0,1),a,b,mu,Sigma,lambda = c(-2,0,1),dist = "SN")
value4 = momentsTMD(c(2,0,1),a,b,mu,Sigma,lambda = c(-2,0,1),tau = 1,dist = "ESN")

\dontrun{
#T cases with kappa scalar (all moments up to 3)
value5 = momentsTMD(3,a,b,mu,Sigma,nu = 7,dist = "t")
value6 = momentsTMD(3,a,b,mu,Sigma,lambda = c(-2,0,1),nu = 7,dist = "ST")
value7 = momentsTMD(3,a,b,mu,Sigma,lambda = c(-2,0,1),tau = 1,nu = 7,dist = "EST")
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