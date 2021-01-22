\name{robustRegBS}
\alias{robustRegBS}
\title{Robust Fitting of Linear Models using Bisquare Psi Function}
\description{
Using iteratively reweighted least squares (IRLS), the function calculates the optimal weights to perform m-estimator or bounded influence regression.  Returns robust beta estimates, mean squared error (MSE) and prints robust ANOVA table.}
\usage{
robustRegBS(formula,data,tune=4.685,m=TRUE,max.it=1000,tol=1e-5,anova.table=FALSE)}
\details{
M-estimates of beta should be used when evaluating least squares estimates of beta and diagnostics show outliers.  Least squares estimates of beta should be used as starting points to achieve convergence.

Bounded influence estimates of beta should be used when evaluating least squares estimates of beta and diagnostics show large values of the "Hat Matrix" diagonals and outliers.}
\arguments{
\item{formula}{Model}
\item{data}{A data frame containing the variables in the model.}
\item{tune}{Tuning Constant.  Default value of 4.685 is 95\% asymptotically efficient against outliers}
\item{m}{If \code{TRUE}, calculates m estimates of beta.  If \code{FALSE}, calculates bounded influence estimates of beta}
\item{max.it}{Maximum number of iterations to achieve convergence in IRLS algorithm}
\item{tol}{Tolerance level in determining convergence}
\item{anova.table}{If \code{TRUE}, prints robust ANOVA table}
}
                                                                                                                                              \author{Ian M. Johnson}
\note{Original package written in 2006}
\references{Tukey,

Birch, Robust F-Test, 1983}
\seealso{robustRegH()}
\examples{
data(stackloss)
robustRegBS(stack.loss~Air.Flow+Water.Temp,data=stackloss)

#If X matrix contained large values of H matrix (high influence points)
robustRegBS(stack.loss~Air.Flow+Water.Temp,data=stackloss,m=FALSE)
}
\keyword{regression}

