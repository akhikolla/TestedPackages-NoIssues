\name{momentsFMD}
\alias{momentsFMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Moments for folded multivariate distributions
}
\description{
It computes the kappa-th order moments for the folded \code{p}-variate Normal, Skew-normal (SN), Extended Skew-normal (ESN) and Student's t-distribution. It also output other lower moments involved in the recurrence approach.
}
\usage{
momentsFMD(kappa,mu,Sigma,lambda = NULL,tau = NULL,nu = NULL,dist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kappa}{
moments vector of length \eqn{p}. All its elements must be integers greater or equal to \eqn{0}. For the Student's-t case, \code{kappa} can be a scalar representing the order of the moment.
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
It represents the degrees of freedom for the Student's t-distribution. Must be an integer greater than 1.
}
 \item{dist}{
represents the folded distribution to be computed. The values are \code{normal}, \code{SN} , \code{ESN} and \code{t}  for the doubly truncated Normal, Skew-normal, Extended Skew-normal and Student's t-distribution respectively.
}
}
\details{
Univariate case is also considered, where \code{Sigma} will be the variance \eqn{\sigma^2}.
}
\value{
A data frame containing \eqn{p+1} columns. The \eqn{p} first containing the set of combinations of exponents summing up to \code{kappa} and the last column containing the the expected value. Normal cases (ESN, SN and normal) return \code{prod(kappa)+1} moments while the Student's t-distribution case returns all moments of order up to \code{kappa}. See example section.
}

\references{
Arellano-Valle, R. B. & Genton, M. G. (2010). Multivariate extended skew-t distributions and related
families. Metron, 68(3), 201-234.

Galarza C.E., Matos L.A., Dey D.K. & Lachos V.H. (2019) On Moments of Folded and Truncated Multivariate Extended Skew-Normal Distributions. Technical report. ID 19-14. University of Connecticut <\url{https://stat.uconn.edu/tech-reports-2019/}>.

Kan, R., & Robotti, C. (2017). On moments of folded and truncated multivariate normal distributions. Journal of Computational and Graphical Statistics, 26(4), 930-934.
}
\author{Christian E. Galarza <\email{cgalarza88@gmail.com}> and
Victor H. Lachos <\email{hlachos@uconn.edu}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}
\note{
Degrees of freedom must be a positive integer. If \code{nu >= 300}, Normal case is considered."
}

\section{Warning }{For the Student-t cases, including ST and EST, \code{kappa}-\eqn{th} order moments exist only for  \code{kappa < nu}.}

\seealso{
\code{\link{meanvarFMD}}, \code{\link{onlymeanTMD}},\code{\link{meanvarTMD}},\code{\link{momentsTMD}}, \code{\link{dmvSN}},\code{\link{pmvSN}},\code{\link{rmvSN}}, \code{\link{dmvESN}},\code{\link{pmvESN}},\code{\link{rmvESN}}, \code{\link{dmvST}},\code{\link{pmvST}},\code{\link{rmvST}}, \code{\link{dmvEST}},\code{\link{pmvEST}},\code{\link{rmvEST}}
}

\examples{
mu = c(0.1,0.2,0.3)
Sigma = matrix(data = c(1,0.2,0.3,0.2,1,0.4,0.3,0.4,1),
               nrow = length(mu),ncol = length(mu),byrow = TRUE)
value1 = momentsFMD(c(2,0,1),mu,Sigma,dist="normal")
value2 = momentsFMD(3,mu,Sigma,dist = "t",nu = 7)
value3 = momentsFMD(c(2,0,1),mu,Sigma,lambda = c(-2,0,1),dist = "SN")
value4 = momentsFMD(c(2,0,1),mu,Sigma,lambda = c(-2,0,1),tau = 1,dist = "ESN")

#T case with kappa vector input
value5 = momentsFMD(c(2,0,1),mu,Sigma,dist = "t",nu = 7)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Folded}% use one of  RShowDoc("KEYWORDS")
\keyword{Multivariate}% __ONLY ONE__ keyword per line
\keyword{Student's t}
\keyword{Normal}
\keyword{Skew}
\keyword{Extended}
