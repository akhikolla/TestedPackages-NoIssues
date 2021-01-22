\name{grpCox}
\alias{grpCox}
\title{Fit a penalized Cox model.}
\description{Fit the regularization paths for Cox models with grouped covariates.}
\usage{
grpCox(X, y, g, m, penalty=c("glasso", "gSCAD", "gMCP"), lambda=NULL, 
nlambda=100, rlambda=NULL, gamma=switch(penalty, gSCAD = 3.7, 3), 
standardize=TRUE, thresh=1e-3, maxit=1e+4)
}

\arguments{
  \item{X}{The design matrix.}
  \item{y}{The response vector includes time corresponding to failure/censor times, and   status indicating failure (1) or censoring (0).}
  \item{g}{A vector indicating the group structure of the covariates. It can be unordered groups.}
  \item{m}{Group multipliers. Default is the square root of group size.}
  \item{penalty}{The penalty to be applied to the model. It is one of \code{glasso}, \code{gSCAD}, or \code{gMCP}.}
  \item{lambda}{A user supplied sequence of \code{lambda} values. If it is left unspecified, and the function automatically computes a grid of lambda values.}
  \item{nlambda}{The number of \code{lambda} values to use in the regularization path.  Default is 100.}
  \item{rlambda}{Smallest value for lambda, as a fraction of the maximum lambda, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero). The default depends on the sample size relative to the number of covariates. If sample size>#covariates, the default is {0.001}, close to zero.  If sample size>#covariates, the default is 0.05.}
  \item{gamma}{Tuning parameter of the group SCAD/MCP penalty.  Default is 3.7 for SCAD and 3 for MCP.}
  \item{standardize}{Logical flag for variable standardization prior to fitting the model.}
  \item{thresh}{Convergence threshold for one-step coordinate descent. Defaults value is 1E-7.}
  \item{maxit}{Maximum number of passes over the data for all lambda values; default is 1E+5.}
}

\details{
The the group SCAD (\code{gSCAD}) and group MCP (\code{gMCP}) formulations have been presented in Wang et. al 2007, Huang et. al 2012.
}
\value{
   \item{aBetaSTD}{A standardized coefficient matrix whose columns correspond to nlambda values of lambda.}
   \item{aBetaO}{A coefficient matrix (without standardization) whose columns correspond to nlambda values of lambda.}
   \item{lambda}{The lambda values used.}
   \item{ll}{The log likelihood values.}
   \item{g}{A vector indicating the group structure of the covariates.}
}

\references{
Wang, L., Chen, G., and Li, H. Group SCAD regression analysis for microarray time course gene expression
data. Bioinformatics 23.12 (2007), pp. 1486-1494.

Huang, J., Breheny, P., and Ma, S. A selective review of group selection in high-dimensional models. Statistical Science 27.4 (2012), pp. 481-499.
}

\author{
Xuan Dang <\email{xuandang11289@gmail.com}>
}

\examples{
set.seed(200)
N <- 50
p <- 9
x <- matrix(rnorm(N * p), nrow = N)
beta <- c(.65,.65,0,0,.65,.65,0,.65,0)
hx <- exp(x \%*\% beta) 
ty <- rexp(N,hx) 
tcens <- 1 - rbinom(n=N, prob = 0.2, size = 1)
y <- data.frame(illt=ty, ills=tcens)
names(y) <- c("time", "status")

g <- c(1,1,2,2,3,3,2,3,2)
m <- c(sqrt(2),sqrt(4),sqrt(3))

fit <- grpCox(x,y,g,m,penalty="glasso")
plot.gCoef(fit$aBetaO, fit$g, fit$lambda)
}

\keyword{Cox models}
\keyword{group regularization}
