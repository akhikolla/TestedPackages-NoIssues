\name{cv.grpCoxOverlap}
\alias{cv.grpCoxOverlap}
\title{
Cross-validation for grpCoxOverlap
}
\description{
Does k-fold cross-validation for grpCoxOverlap
}
\usage{
cv.grpCoxOverlap(X0, y, group, penalty=c("glasso", "gSCAD", "gMCP"), 
lambda=NULL, nlambda=100, rlambda=NULL,gamma=switch(penalty, SCAD = 3.7, 3), 
standardize=TRUE, thresh=1e-3, maxit=1e+4, nfolds=10, foldid=NULL, 
returnLatent=TRUE)
}

\arguments{
  \item{X0}{The design matrix.}
  \item{y}{The response vector includes time corresponding to failure/censor times, and   status indicating failure (1) or censoring (0).}
  \item{group}{A list of groups, each includes indices of covariates in the group.}
  \item{penalty}{The penalty to be applied to the model. It is one of \code{glasso}, \code{gSCAD}, or \code{gMCP}.}
  \item{lambda}{A user supplied sequence of \code{lambda} values. If it is left unspecified, and the function automatically computes a grid of lambda values.}
  \item{nlambda}{The number of \code{lambda} values to use in the regularization path.  Default is 100.}
  \item{rlambda}{Smallest value for lambda, as a fraction of the maximum lambda, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero). The default depends on the sample size relative to the number of covariates. If sample size>#covariates, the default is {0.001}, close to zero.  If sample size>#covariates, the default is 0.05.}
  \item{gamma}{Tuning parameter of the group SCAD/MCP penalty. Default is 3.7 for SCAD and 3 for MCP.}
  \item{standardize}{Logical flag for variable standardization prior to fitting the model.}
  \item{thresh}{Convergence threshold for one-step coordinate descent. Defaults value is 1E-7.}
  \item{maxit}{Maximum number of passes over the data for all lambda values; default is 1E+5.}
  \item{nfolds}{The number of cross-validation folds. Default is 10.}
  \item{foldid}{An optional vector of values between 1 and nfolds identifying what fold each observation is in.}
  \item{returnLatent}{Return the coefficient matrix in latent space. Default is TRUE.}
}

               
\value{
    \item{aBetaLatent}{A coefficient matrix whose columns correspond to nlambda values of lambda in latent space.}
   \item{aBetaOri}{A coefficient matrix whose columns correspond to nlambda values of lambda in original space.}
  \item{mBetaLatent}{The coefficient in latent space gives maximum log-likelihood value using the first cross-validation method.}
  \item{mBetaOri}{The coefficient in original space gives maximum log-likelihood value using the first cross-validation method.}
  \item{pBetaLatent}{The coefficient in latent space gives maximum log-likelihood value using the penalized cross-validation method.}
  \item{pBetaOri}{The coefficient in original space gives maximum log-likelihood value using the penalized cross-validation method.}
  \item{fit}{A matrix includes lambda value, the mean cross-validation error.}
  \item{lambda}{The lambda values used.}
  \item{group}{A list of groups, each includes indices of covariates in the group.}
  \item{glatent}{A vector indicating the group structure of the covariates in latent space.}
  \item{cvmax}{The maximum value of log likelihood.}
  \item{lambda.max}{The value of lambda corresponds to the maximum value of log likelihood using the first cross-validation method.}
  \item{lambda.pcvl}{The value of lambda corresponds to the maximum value of log likelihood using the penalized cross-validation method.}
}

\references{
Verweij PJ, Houwelingen HC. Cross-validation in survival analysis. Statistics in Medicine 1993; 12(24): 385-395.

Ternes N, Rotolo F, Michiels S. Empirical extensions of the lasso penalty to reduce the false discovery rate in highdimensional Cox regression models. Statistics in Medicine 2016; 35(15): 2561-73.
}

\author{
Xuan Dang <\email{xuandang11289@gmail.com}>
}

\examples{
set.seed(100001)
N <- 50
p <- 6
times <- 1:p
rho <- 0.5
H <- abs(outer(times, times, "-"))
C <- 1 * rho^H
C[cbind(1:p, 1:p)] <- C[cbind(1:p, 1:p)] 
sigma <- matrix(C,p,p)
mu <- rep(0,p)
x <- mvrnorm(n=N, mu, sigma)

beta <- c(0, .8, 1, 2, 1, 0)
hx <- exp(x \%*\% beta) 
ty <- rexp(N,hx) 
tcens <- 1 - rbinom(n=N, prob = 0.2, size = 1)
y <- data.frame(illt=ty, ills=tcens)
names(y) <- c("time", "status")

group <- list(g1 = c(1,2,3,4), g2 = c(1,2,6), g3 = c(2,3), 
              g4 = c(4,5), g5 = c(5))
cvfit <- cv.grpCoxOverlap(x, y, group, penalty="glasso", nlambda=50)
plot.llCV(cvfit)
}

\keyword{Cox models}
\keyword{group regularization}
\keyword{overlapping group}