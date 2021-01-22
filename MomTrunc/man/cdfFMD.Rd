\name{cdfFMD}
\alias{cdfFMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cumulative distribution function for folded multivariate distributions
}
\description{
It computes the cumulative distribution function on \code{x} for a folded \code{p}-variate Normal, Skew-normal (SN), Extended Skew-normal (ESN) and Student's t-distribution.
}
\usage{
cdfFMD(x,mu,Sigma,lambda = NULL,tau = NULL,dist,nu = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of length \eqn{p} where the cdf is evaluated.
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
  \item{dist}{
represents the folded distribution to be computed. The values are \code{normal}, \code{SN} , \code{ESN} and \code{t}  for the doubly truncated Normal, Skew-normal, Extended Skew-normal and Student's t-distribution respectively.
}
  \item{nu}{
It represents the degrees of freedom for the Student's t-distribution.
}
}
\details{
Normal case by default, i.e., when \code{dist} is not provided. Univariate case is also considered, where \code{Sigma} will be the variance \eqn{\sigma^2}.
}
\value{
It returns the distribution value for a single point \code{x}.
}
\references{
Chakraborty, A. K., & Chatterjee, M. (2013). On multivariate folded normal distribution. Sankhya B, 75(1), 1-15.

	Galarza C.E., Matos L.A., Dey D.K. & Lachos V.H. (2019) On Moments of Folded and Truncated Multivariate Extended Skew-Normal Distributions. Technical report. ID 19-14. University of Connecticut <\url{https://stat.uconn.edu/tech-reports-2019/}>.
	
Kan R. & Robotti C. (2017) On Moments of Folded and Truncated Multivariate Normal Distributions, Journal of Computational and Graphical Statistics, 26:4, 930-934.
}

\author{Christian E. Galarza <\email{cgalarza88@gmail.com}> and
Victor H. Lachos <\email{hlachos@uconn.edu}>

Maintainer: Christian E. Galarza <\email{cgalarza88@gmail.com}>
}
\note{
Degrees of freedom must be a positive integer. If \code{nu >= 200}, Normal case is considered."
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{momentsFMD}}, \code{\link{meanvarFMD}}
}
\examples{
mu = c(0.1,0.2,0.3,0.4)
Sigma = matrix(data = c(1,0.2,0.3,0.1,0.2,1,0.4,-0.1,0.3,0.4,1,0.2,0.1,-0.1,0.2,1),
               nrow = length(mu),ncol = length(mu),byrow = TRUE)
cdfFMD(x = c(0.5,0.2,1.0,1.3),mu,Sigma,dist="normal")
cdfFMD(x = c(0.5,0.2,1.0,1.3),mu,Sigma,dist = "t",nu = 4)
cdfFMD(x = c(0.5,0.2,1.0,1.3),mu,Sigma,lambda = c(-2,0,2,1),dist = "SN")
cdfFMD(x = c(0.5,0.2,1.0,1.3),mu,Sigma,lambda = c(-2,0,2,1),tau = 1,dist = "ESN")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Folded}% use one of  RShowDoc("KEYWORDS")
\keyword{Multivariate}% __ONLY ONE__ keyword per line
\keyword{Normal}
\keyword{Skew}
\keyword{Extended}
\keyword{Student's t}
