\name{SIRladle}
\alias{SIRladle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ladle Estimate for SIR
}
\description{
In the supervised dimension reduction context with response y and explaining variables x, this functions provides the ladle estimate 
for the dimension of the central subspace for SIR.
}
\usage{
SIRladle(X, y, h = 10, n.boot = 200, 
         ncomp = ifelse(ncol(X) > 10, floor(ncol(X)/log(ncol(X))), ncol(X) - 1), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{numeric data matrix.}
  \item{y}{numeric response vector.}
  \item{h}{number of slices in SIR.}
  \item{n.boot}{number of bootstrapping samples to be used.}
  \item{ncomp}{The number of components among which the ladle estimator is to be searched. The default here follows
  the recommendation of Luo and Li 2016.}
  \item{\dots}{arguments passed on to \code{\link{quantile}}.}
}
\details{
The idea here is that the eigenvalues of the SIR-M matrix are of the form  \eqn{\lambda_1 \geq ... \geq \lambda_k > 0 =  ... = 0}{lambda_1 >= ... >= lambda_k > 0 =  ... = 0} and the eigenvectors
of the non-zero eigenvalue span the central subspace.

The ladle estimate for k for this purpose combines the values of the 
scaled eigenvalues and the variation of the eigenvectors based on bootstrapping. The idea there is that for distinct eigenvales the variation of the eigenvectors
is small and for equal eigenvalues the corresponding eigenvectors have large variation.

This measure is then computed assuming k=0,..., \code{ncomp} and the ladle estimate for k is the value where the measure takes its minimum. 

}
\value{
A list of class ladle containing:
  \item{method}{the string SIR.}
  \item{k}{the estimated value of k.}
  \item{fn}{vector giving the measures of variation of the eigenvectors using the bootstrapped eigenvectors for the different number of components.}
  \item{phin}{normalized eigenvalues of the M matrix in the SIR case.}
  \item{gn}{the main criterion for the ladle estimate - the sum of fn and phin. k is the value where gn takes its minimum}
  \item{lambda}{the eigenvalues of the M matrix in the SIR case.}
  \item{W}{the transformation matrix to supervised components.}
  \item{S}{data matrix with the centered supervised components.}
  \item{MU}{the location of the data which was substracted before calculating the supervised components.}
  \item{data.name}{the name of the data for which the ladle estimate was computed.}
}
\references{
\cite{Luo, W. and Li, B. (2016), Combining Eigenvalues and Variation of Eigenvectors for Order Determination, Biometrika, 103. 875--887. <doi:10.1093/biomet/asw051>}
}
\author{
Klaus Nordhausen
}
\seealso{
\code{\link{ladleplot}}
}
\examples{
n <- 1000
X <- cbind(rnorm(n), rnorm(n), rnorm(n), rnorm(n), rnorm(n))
eps <- rnorm(n, sd=0.02)
y <- 4*X[,1] + 2*X[,2] + eps


test <- SIRladle(X, y)
test
summary(test)
plot(test)
pairs(cbind(y, components(test)))
ladleplot(test)
ladleplot(test, crit = "fn")
ladleplot(test, crit = "lambda")
ladleplot(test, crit = "phin")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")
