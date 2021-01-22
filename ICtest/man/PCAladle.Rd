\name{PCAladle}
\alias{PCAladle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ladle Estimate for PCA
}
\description{
For p-variate data, the Ladle estimate for PCA assumes that the last p-k eigenvalues are equal. Combining information from the eigenvalues and eigenvectors
of the covariance matrix the ladle estimator yields an estimate for k. 
}
\usage{
PCAladle(X, n.boot = 200, 
         ncomp = ifelse(ncol(X) > 10, floor(ncol(X)/log(ncol(X))), ncol(X) - 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{numeric data matrix.}
  \item{n.boot}{number of bootstrapping samples to be used.}
  \item{ncomp}{The number of components among which the ladle estimator is to be searched. The default here follows
  the recommendation of Luo and Li 2016.}
}
\details{
The model here assumes that the eigenvalues of the covariance matrix are of the form \eqn{\lambda_1 \geq ... \geq \lambda_{k} > \lambda_{k+1} =  ... = \lambda_p}{lambda_1 >= ... >= lambda_k > lambda_k+1 =  ... = \lambda_p}
and the goal is to estimate the value of k. The ladle estimate for this purpose combines the values of the 
scaled eigenvalues and the variation of the eigenvectors based on bootstrapping. The idea there is that for distinct eigenvales the variation of the eigenvectors
is small and for equal eigenvalues the corresponding eigenvectors have large variation.

This measure is then computed assuming k=0,..., \code{ncomp} and the ladle estimate for k is the value where the measure takes its minimum. 

}

\value{
A list of class ladle containing:
  \item{method}{the string PCA.}
  \item{k}{the estimated value of k.}
  \item{fn}{vector giving the measures of variation of the eigenvectors using the bootstrapped eigenvectors for the different number of components.}
  \item{phin}{normalized eigenvalues of the covariance matrix.}
  \item{gn}{the main criterion for the ladle estimate - the sum of fn and phin. k is the value where gn takes its minimum}
  \item{lambda}{the eigenvalues of the covariance matrix.}
  \item{W}{the transformation matrix to the principal components.}
  \item{S}{data matrix with the centered principal components.}
  \item{MU}{the location of the data which was substracted before calculating the principal components.}
  \item{data.name}{the name of the data for which the ladle estimate was computed.}
}
\references{
\cite{Luo, W. and Li, B. (2016), Combining Eigenvalues and Variation of Eigenvectors for Order Determination, Biometrika, 103, 875--887. <doi:10.1093/biomet/asw051>}
}
\author{
Klaus Nordhausen
}
\seealso{
\code{\link{ladleplot}}
}
\examples{
n <- 1000
Y <- cbind(rnorm(n, sd=2), rnorm(n,sd=2), rnorm(n), rnorm(n), rnorm(n), rnorm(n))

testPCA <- PCAladle(Y)
testPCA
summary(testPCA)
plot(testPCA)
ladleplot(testPCA)
ladleplot(testPCA, crit = "fn")
ladleplot(testPCA, crit = "lambda")
ladleplot(testPCA, crit = "phin")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")
