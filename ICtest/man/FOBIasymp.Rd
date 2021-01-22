\name{FOBIasymp}
\alias{FOBIasymp}

\title{
Testing for the Number of Gaussian Components in NGCA or ICA Using FOBI 
}
\description{
In non-gaussian component analysis (NGCA) and independent components analysis (ICA) gaussian components are considered as uninteresting.
The function tests, based on FOBI, if there are \code{p-k} gaussian components where \code{p} is the dimension of the data.
The function offers three different test versions.
}
\usage{
FOBIasymp(X, k, type = "S3", model = "NGCA", method = "satterthwaite")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{numeric data matrix.}
  \item{k}{the number of non-gaussian components under the null.}
  \item{type}{which of the three tests to perform. Options are \code{"S1"}, \code{"S2"} and \code{"S3"}. For the differences see the details section.}
  \item{model}{What is the underlying assumption of the non-gaussian parts. Options are general \code{"NGCA"} model and \code{"ICA"} model.}
  \item{method}{if \code{type = "S1"} the teststatistic has as limiting distribution a weighted sum of chisquare distributions. To compute the p-value then the 
  function used is \code{\link[survey]{pchisqsum}}. The \code{method} argument specifies which method  \code{\link[survey]{pchisqsum}} uses for the computation.
  Options are \code{"satterthwaite"}, \code{"integration"} and  \code{"saddlepoint"}.}
}

\details{
The function jointly diagonalizes the regular covariance and the matrix of fourth moments. Note however that in this case the matrix of fourth moments
is not made consistent under the normal model by dividing it by \eqn{p+2}, as for example done by the function \code{\link[ICS]{cov4}} where \eqn{p} denotes the dimension
of the data. Therefore 
the eigenvalues of this generalized eigenvector-eigenvalue problem which correspond to normally distributed components should be \code{p+2}.

Given eigenvalues \eqn{d_1,...,d_p} the function thus orders the components in decending order according to the values of \eqn{(d_i-(p+2))^2}.

Under the null it is then assumed that the first \code{k} interesting components are mutually independent and non-normal and the last \code{p-k} are gaussian.

Three possible tests are then available to test this null hypothesis for a sample of size n:

\enumerate{
\item \code{type="S1"}: The test statistic T is the variance of the last p-k eigenvalues around p+2:
\deqn{T = n \sum_{i=k+1}^p (d_i-(p+2))^2}{T = n sum_{i=k+1}^p (d_i-(p+2))^2,}
the limiting distribution of which under the null is the sum of two weighted chisquare distributions with weights: 

\eqn{w_1 = 2 \sigma_1 / (p-k)}{w_1 = 2*sigma_1 / (p-k)} and \eqn{w_2 = 2 \sigma_1 / (p-k)  +  \sigma_2}{w_2 = 2*sigma_1 / (p-k) + sigma_2}.

and degrees of freedom:

\eqn{df_1 = (p-k-1)(p-k+2)/2} and \eqn{df_2 = 1}.

\item \code{type="S2"}: Another possible version for the test statistic is a scaled sum of the variance of the eigenvalues around the mean plus the variance around
the expected value under normality (p+2). Denote \eqn{VAR_{dpk}}{VAR_dpk} as the variance of the last p-k eigenvalues  and  \eqn{VAR2_{dpk}}{VAR2_dpk} as the variance of these eigenvalues around \eqn{p+2}.
Then the test statistic is:
\deqn{T = (n (p-k) VAR_{dpk}) / (2 \sigma_1) + (n VAR2_{dpk}) / (2 \sigma_1 / (p-k) + \sigma_2)}{T = (n (p-k) VAR_dpk) / (2 sigma_1) + (n VAR2_dpk) / (2 sigma_1 / (p-k) + sigma_2).}

This test statistic has a limiting chisquare distribution with \eqn{(p-k-1)(p-q+2)/2 + 1} degrees of freedom.

\item \code{type="S3"}: The third possible test statistic just checks the equality of the last p-k eigenvalues using only the first part of the test statistic of \code{type="S2"}.
The test statistic is then:
\deqn{T = (n (p-k) VAR_{dpk}) / (2 \sigma_1)}{T = (n (p-k) VAR_dpk) / (2 sigma_1)}

and has a limiting chisquare distribution with \eqn{(p-k-1)(p-q+2)/2} degrees of freedom.
}

The constants \eqn{\sigma_1}{sigma_1} and \eqn{\sigma_2}{sigma_2} depend on the underlying model assumptions as specified in argument \code{model} and are estimated from the data.


}
\value{
A list of class ictest inheriting from class htest containing:
  \item{statistic}{the value of the test statistic.}
  \item{p.value }{the p-value of the test.}
  \item{parameter }{the degrees of freedom of the test or the degrees of freedoms and the corresponding weights of the test in case the test
  has as its limiting distribution a weighted sum of chisquare distributions.}
  \item{method}{character string denoting which test was performed.}
  \item{data.name}{character string giving the name of the data.}
  \item{alternative}{character string specifying the alternative hypothesis.}
  \item{k}{the number or non-gaussian components used in the testing problem.}
  \item{W}{the transformation matrix to the independent components. Also known as unmixing matrix.}
  \item{S}{data matrix with the centered independent components.}
  \item{D}{the underlying FOBI eigenvalues.}
  \item{MU}{the location of the data which was substracted before calculating the independent components.}
  \item{sigma1}{the asymptotic constant sigma1 needed for the asymptotic test(s).}
  \item{sigma2}{the asymptotic constant sigma2 needed for the asymptotic test(s).}
  \item{type}{the value of \code{type}.}
  \item{model}{the value of \code{model}.}
}
\references{
\cite{Nordhausen, K., Oja, H. and Tyler, D.E. (2016), Asymptotic and bootstrap tests for subspace dimension, <https://arxiv.org/abs/1611.04908>.}

\cite{Nordhausen, K., Oja, H., Tyler, D.E. and Virta, J. (2017), Asymptotic and Bootstrap Tests for the Dimension of the Non-Gaussian Subspace, Signal Processing Letters, 24, 887--891. <doi:10.1109/LSP.2017.2696880 >.}
}
\author{
Klaus Nordhausen
}


\seealso{
\code{\link[JADE]{FOBI}},  \code{\link{FOBIboot}}
}
\examples{
n <- 1500
S <- cbind(runif(n), rchisq(n, 2), rexp(n), rnorm(n), rnorm(n), rnorm(n))
A <- matrix(rnorm(36), ncol = 6)
X <- S \%*\% t(A)

FOBIasymp(X, k = 2)
FOBIasymp(X, k = 3, type = "S1")
FOBIasymp(X, k = 0, type = "S2", model = "ICA")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
\keyword{ multivariate }% __ONLY ONE__ keyword per line
