\name{ladle}
\alias{ladle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ladle estimate for an arbitrary matrix
}
\description{
The ladle estimates the rank of a symmetric matrix \eqn{S} by combining the classical screeplot with an estimate of the rank from the bootstrap eigenvector variability of \eqn{S}.
}
\usage{
ladle(x, S, n.boots = 200, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{n} x \code{p} data matrix.}
  \item{S}{Function for computing a \code{q} x \code{q} symmetric matrix from the data \code{x}.}
  \item{n.boots}{The number of bootstrap samples.}
  \item{...}{Furhter parameters passed to \code{S}}
}
\details{
Assume that the eigenvalues of the population version of \code{S} are \eqn{\lambda_1 >= ... >= \lambda_k > \lambda_k+1 = ... = \lambda_p}{lambda_1 >= ... >= lambda_k > lambda_k+1 = ... = lambda_p}. The ladle estimates the true value of \eqn{k} (for example the rank of \code{S}) by combining the classical screeplot with estimate of \eqn{k} from the bootstrap eigenvector variability of \code{S}.

For applying the ladle to either PCA, FOBI or SIR, see the dedicated functions \code{\link{PCAladle}}, \code{\link{FOBIladle}}, \code{\link{SIRladle}}.
}
\value{
A list of class \code{ladle} containing:
\item{method }{The string ``general''.}
\item{k }{The estimated value of k.}
\item{fn }{A vector giving the measures of variation of the eigenvectors using the bootstrapped eigenvectors for the different number of components.}
\item{phin }{The normalized eigenvalues of the S matrix.}
\item{gn }{The main criterion for the ladle estimate - the sum of fn and phin. k is the value where gn takes its minimum.}
\item{lambda}{The eigenvalues of the covariance matrix.}
\item{data.name}{The name of the data for which the ladle estimate was computed.}
}
\references{
Luo, W. and Li, B. (2016), Combining Eigenvalues and Variation of Eigenvectors for Order Determination, Biometrika, 103. 875-887. <doi:10.1093/biomet/asw051>
}
\author{
Joni Virta
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{PCAladle}}, \code{\link{FOBIladle}}, \code{\link{SIRladle}}
}
\examples{
# Function for computing the left CCA matrix
S_CCA <- function(x, dim){
  x1 <- x[, 1:dim]
  x2 <- x[, -(1:dim)]
  stand <- function(x){
    x <- as.matrix(x)
    x <- sweep(x, 2, colMeans(x), "-")
    eigcov <- eigen(cov(x), symmetric = TRUE)
    x\%*\%(eigcov$vectors\%*\%diag((eigcov$values)^(-1/2))\%*\%t(eigcov$vectors))
  }
  
  x1stand <- stand(x1)
  x2stand <- stand(x2)
  
  crosscov <- cov(x1stand, x2stand)
  
  tcrossprod(crosscov)
}

# Toy data with two canonical components
n <- 200
x1 <- matrix(rnorm(n*5), n, 5)
x2 <- cbind(x1[, 1] + rnorm(n, sd = sqrt(0.5)),
            -1*x1[, 1] + x1[, 2] + rnorm(n, sd = sqrt(0.5)),
            matrix(rnorm(n*3), n, 3))
x <- cbind(x1, x2)

# The ladle estimate
ladle_1 <- ladle(x, S_CCA, dim = 5)
ladleplot(ladle_1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")

