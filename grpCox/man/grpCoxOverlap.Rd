\name{grpCoxOverlap}
\alias{grpCoxOverlap}
\title{Fit a penalized regression path with overlapping grouped covariates.}
\description{Fit the regularization paths for Cox's models with overlapping grouped covariates.}
\usage{
grpCoxOverlap(X0, y, group, penalty=c("glasso", "gSCAD", "gMCP"), 
lambda=NULL, nlambda=100, rlambda=NULL, gamma=switch(penalty, gSCAD = 3.7, 3),
standardize = TRUE, thresh=1e-3, maxit=1e+4, returnLatent=TRUE)
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
  \item{returnLatent}{Return the coefficient matrix in latent space. Default is TRUE.}
}

\details{
The the group SCAD (\code{gSCAD}) and group MCP (\code{gMCP}) formulations have been presented in Wang et. al 2007, Huang et. al 2012.

The method based on the latent group approach (Jacob et al. 2009, Obozinski et al. 2011.)
}

\value{
   \item{aBetaLatent}{A coefficient matrix whose columns correspond to nlambda values of lambda in latent space.}
   \item{aBetaOri}{A coefficient matrix whose columns correspond to nlambda values of lambda in original space.}
   \item{lambda}{The lambda values used.}
   \item{ll}{The log likelihood values.}
   \item{group}{A list of groups, each includes indices of covariates in the group.}
   \item{glatent}{A vector indicating the group structure of the covariates in latent space.}
}

\references{
Wang, L., Chen, G., and Li, H. Group SCAD regression analysis for microarray time course gene expression
data. Bioinformatics 23.12 (2007), pp. 1486-1494.

Huang, J., Breheny, P., and Ma, S. A selective review of group selection in high-dimensional models." Statistical Science 27.4 (2012), pp. 481-499.

Jacob, L., Obozinski, G., and Vert, J. P. (2009, June). Group lasso with overlap and graph lasso. In Proceedings of the 26th annual international conference on machine learning, ACM: 433-440.

Obozinski, G., Jacob, L., and Vert, J. P. (2011). Group lasso with overlaps: the latent group lasso approach.
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

group <- list(g1 = c(1,2,3,4), g2 = c(1,2,6), g3 = c(2,3), g4 = c(4,5), g5 = c(5))
fit <- grpCoxOverlap(x, y, group, penalty="glasso", nlambda=50)
# plot the coefficient values in latent space
plot.gCoef(fit$aBetaLatent, fit$glatent, fit$lambda)
# plot the coefficient values in original space
plot.Coef(fit$aBetaOri, fit$lambda)
}

\keyword{Cox models}
\keyword{group regularization}
\keyword{overlapping group}