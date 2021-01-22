\name{MomTrunc-package}
\alias{MomTrunc-package}
\alias{MomTrunc}
\docType{package}
\title{
\packageTitle{MomTrunc}
}
\description{
\packageDescription{MomTrunc}
}
\details{
Probabilities can be computed using the functions \code{\link{pmvSN}} and \code{\link{pmvESN}} for the normal cases SN and ESN and, \code{\link{pmvST}} and \code{\link{pmvEST}} for the t cases ST and EST respectively, which offer the option to return the logarithm in base 2 of the probability, useful when the true probability is too small for the machine precision. These functions above use methods in Genz (1992) through the \code{mvtnorm} package (linked direclty to our \code{C++} functions) and Cao et.al. (2019) through the package \code{tlrmvnmvt}. For the double truncated Student-t cases SUT, EST, ST and T, decimal degrees of freedom are supported. Computation of arbitrary moments are based in the works of Galarza et.al. (2019) and Kan & Robotti (2017). Reference for the family of selection-elliptical distributions in this package can be found in Arellano-Valle & Genton (2005).
}
\author{
\packageAuthor{MomTrunc}

Maintainer: \packageMaintainer{MomTrunc}
}
\references{
Arellano-Valle, R. B. & Genton, M. G. (2005). On fundamental skew distributions. Journal of Multivariate
Analysis, 96, 93-116.

Cao, J., Genton, M. G., Keyes, D. E., & Turkiyyah, G. M. (2019) "Exploiting Low
	Rank Covariance Structures for Computing High-Dimensional Normal and Student-
	t Probabilities" <\url{https://marcgenton.github.io/2019.CGKT.manuscript.pdf}>.
	
Galarza C.E., Matos L.A., Dey D.K. & Lachos V.H. (2019) On Moments of Folded and Truncated Multivariate Extended Skew-Normal Distributions. Technical report. ID 19-14. University of Connecticut <\url{https://stat.uconn.edu/tech-reports-2019/}>.

Genz, A., "Numerical computation of multivariate normal probabilities," Journal of Computational and Graphical Statistics, 1, 141-149 (1992) <doi:10.1080/10618600.1992.10477010>.

Kan, R., & Robotti, C. (2017). On moments of folded and truncated multivariate normal distributions. Journal of Computational and Graphical Statistics, 26(4), 930-934.
}

\keyword{Truncated}% use one of  RShowDoc("KEYWORDS")
\keyword{Multivariate}% __ONLY ONE__ keyword per line
\keyword{Student's t}
\keyword{Normal}
\keyword{Skew}
\keyword{Extended}
\keyword{Unified}
\keyword{Selection}

\seealso{
\code{\link{onlymeanTMD}},\code{\link{meanvarTMD}},\code{\link{momentsTMD}},\code{\link{dmvSN}},\code{\link{pmvSN}},\code{\link{rmvSN}},\code{\link{dmvST}},\code{\link{pmvST}},\code{\link{rmvST}}
}
\examples{
a = c(-0.8,-0.7,-0.6)
b = c(0.5,0.6,0.7)
mu = c(0.1,0.2,0.3)
Sigma = matrix(data = c(1,0.2,0.3,0.2,1,0.4,0.3,0.4,1),
               nrow = length(mu),ncol = length(mu),byrow = TRUE)

meanvarTMD(a,b,mu,Sigma,dist="normal") #normal case
meanvarTMD(mu = mu,Sigma = Sigma,lambda = c(-2,0,1),dist="SN") #skew normal with NO truncation
meanvarTMD(a,b,mu,Sigma,lambda = c(-2,0,1),nu = 4.87,dist = "ST") #skew t
momentsTMD(3,a,b,mu,Sigma,nu = 4,dist = "t") #t case, all moments or order <=3
}
