\name{FOBIladle}
\alias{FOBIladle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ladle Estimate to Estimate the Number of Gaussian Components in ICA or NGCA
}
\description{
The ladle estimator uses the eigenvalues and eigenvectors of FOBI to estimate the number of Gaussian components in ICA or NGCA.
}
\usage{
FOBIladle(X, n.boot = 200, 
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
The model here assumes that in ICA or NGCA there are k non-gaussian components and p-k gaussian components.
The idea is then to decide which eigenvalues differ from p+2. The ladle estimate for this purpose combines the values of the 
scaled eigenvalues and the variation of the eigenvectors based on bootstrapping. The idea there is that for distinct eigenvales the variation of the eigenvectors
is small and for equal eigenvalues the corresponding eigenvectors have large variation.

This measure is then computed assuming k=0,..., \code{ncomp} and the ladle estimate for k is the value where the measure takes its minimum. 

 
}
\value{
A list of class ladle containing:
  \item{method}{the string FOBI.}
  \item{k}{the estimated number of non-gaussian components.}
  \item{fn}{vector giving the measures of variation of the eigenvectors using the bootstrapped eigenvectors for the different number of components.}
  \item{phin}{normalized eigenvalues of the FOBI matrix.}
  \item{gn}{the main criterion for the ladle estimate - the sum of fn and phin. k is the value where gn takes its minimum}
  \item{lambda}{the eigenvalues of the FOBI matrix.}
  \item{W}{the transformation matrix to the independent components. Also known as unmixing matrix.}
  \item{S}{data matrix with the centered independent components.}
  \item{MU}{the location of the data which was substracted before calculating the independent components.}
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
X <- cbind(rexp(n), rt(n,5), rnorm(n), rnorm(n), rnorm(n), rnorm(n))

test <- FOBIladle(X)
test
summary(test)
plot(test)
ladleplot(test)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }% use one of  RShowDoc("KEYWORDS")
