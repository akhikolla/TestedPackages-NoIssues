\name{repolr-package}
\alias{repolr-package}
\alias{repolr} 
\alias{print.repolr}
\alias{coef.repolr}
\alias{fitted.repolr}
\alias{model.matrix.repolr}
\alias{terms.repolr}
\alias{residuals.repolr}
\alias{smat}
\alias{cmat}
\alias{alphpow}
\alias{hgmat}
\alias{ordgee}
\alias{upalpha}
\alias{icormat}
\alias{potest}
\docType{package}
\title{
Repeated Measures Proportional Odds Logistic Regression using GEE}
\description{
The package allows regression models to be fitted to repeated ordinal scores, for the proportional odds model, using a generalized estimating equation (GEE) methodology. The algorithm estimates the correlation parameter by minimizing the generalized variance of the regression parameters at each step of the fitting algorithm. Parameter estimation is available for the uniform and first-order autoregressive correlation models, for data potentially recorded at irregularly spaced time intervals. A test for the proportional odds assumption is available, as is an option to fit a polynomial model to the the cut-point parameters.}
\usage{
repolr(formula, subjects, data, times, categories, corr.mod = "independence",
         alpha = 0.5, po.test = FALSE, fixed = FALSE, poly = NULL, 
         space = NULL, diffmeth = "analytic", fit.opt = rep(NA, 5))
}
\arguments{
  \item{formula}{a formula, as for other regression models.}
  \item{subjects}{a character string specifying the name of the subject variable.}
  \item{data}{a dataframe in which to interpret the variables occurring in the formula.}
  \item{times}{a vector of times which occur within subject clusters; e.g. for four evenly spaced times c(1, 2, 3, 4).}
  \item{categories}{a numeric variable indicating the number of ordinal score categories.}
  \item{corr.mod}{a character string specifying the correlation structure. The following are permitted: \dQuote{\code{ar1}}, \dQuote{\code{uniform}} and \dQuote{\code{independence}}.}
  \item{alpha}{an initial value for the correlation parameter.}  
  \item{po.test}{a logical variable; if true a score test for proportional odds is reported.}
  \item{fixed}{a logical variable; if true the correlation is fixed at the initial value (\code{alpha}) during model fitting.}
  \item{poly}{a numeric variable indicating the order of the polynomial contrasts used for the cut-point model.}
  \item{space}{a vector indicating the category spacing when fitting the polynomial model; can generally be set to 1:{\code{categories}}}
  \item{diffmeth}{a character string specifying the method used for estimation of \code{alpha}. The folowing are available:\dQuote{\code{analytic}} and \dQuote{\code{numeric}}}
  \item{fit.opt}{a vector of options to control the behaviour of the fitting algorithm.}
}
\details{
The \code{repolr} function fits models for repeated ordinal scores using GEE methodology.

The user is required to specify, as a minimum: (i) a data set name (\code{data}), (ii) a model formula (\code{formula}), (iii) a cluster identification variable (\code{subjects}), (iv) a time variable (\code{time}) and (v) the number of categories used for the response variable (\code{categories}). 

The data may contain records with missing data for either the response variable or the explanatory variables. The response variable must have at least three ordered categories (\emph{K} greater than or equal to 3) and, as \emph{K}-1 cut-point parameters are estimated, an additional intercept term can not be explicitly included in the model formula. A subject variable, which takes integer values from 1 to \emph{N} with no missing values allowed, indicates the data clusters (patients or experimental units) and a time variable indicates the within cluster ordering; times must be ordered integers starting from one and spaced to indicate the relative distance between successive times. For instance, four observations at equally spaced times would be entered as 1, 2, 3 and 4, whereas if the first two observations were separated by half the time interval of the other observations then coding would be 1, 2, 4 and 6. The data must be pre-sorted by time clusters within each subject, and complete, i.e. where data is missing for a particular time it must be indicated as such. the datasets provided with this package provide exemplars of the required data formatting; e.g. \code{\link{HHSpain}} and \code{\link{QoL}}.

The available options for the correlation model (\code{corstr}) are AR1, uniform, fixed and independence, with default setting independence. 

Additionally there are a number of other algorithm related options.

The algorithm is generally robust to the initial value for \code{alpha} (default setting = 0.5), where estimation is rerquired, however a starting value for \code{alpha} can be set. If required the correlation parameter, set via \code{alpha}, can be fixed throughout model fitting, and not updated, by setting the option \code{fixed} to \code{TRUE}.

The partial derivatives of the log of the determinant of the robust variance matrix (generalized variance), with respect to \code{alpha}, can either be determined analytically (\dQuote{\code{analytic}} setting for \code{diffmeth}) or numerically by finite differencing (\dQuote{\code{numeric}} setting for \code{diffmeth}). The latter method is often quicker for complex regression models, or if \emph{K} is large.

Function \code{poly}, specifies the order of fitted orthogonal polynomial contrasts for the cut-point parameters; the default setting for \pkg{repolr} is to fit the complete set of cut-point parameters. Fitting polynomial contrasts can be particularly useful for long ordinal scores, where \emph{K} is large, or where a particular form for the cut-point parameters is preferred; e.g. increasing  uniformly at the extremes of the score scale. The order of the polynomial must be an integer less than  \emph{K}-1. The function requires one additional argument, \code{space}, that indicates the spacing between categories. This is normally set to 1:{\code{categories}} when all categories are observed, but can be modified to provide more realsitic models if for instance some categories are not observed.

The proportional odds assumption can be tested, using a score test, by setting the option \code{po.test} to TRUE.

The \code{fit.opt} argument provides control of the fitting algorithm; the defaults are

\code{c(cmaxit = 10, omaxit = 5, ctol = 0.001, otol = 0.00001, h = 0.01)}.

These control the maximum number of iterations for updating estimates of \code{alpha}, the maximum number of iterations for updating the regression parameters within each of the updating steps for \code{alpha}, the convergence tolerances for estimation of \code{alpha} and the regression parameters, and the interval \code{h} for finite differencing, if the \dQuote{\code{numeric}} option is selected.

Model fitting is implemented via a suite of of functions developed using \pkg{RcppArmadillo} (\pkg{Rcpp}) that construct correlations between derived binary variables at each time-point (\code{smat}) and between time-points (\code{cmat}). Complete (sparse) model covariance matrices are constructed using \code{hgmat} and \code{alphpow}, with model parameter estimation implemented in \code{ordgee}, using current estimates of the inverse of the correlation matrix from \code{icormat}. Function \code{upalpha} provides updates of  correlation parameter estimates, and \code{potest} implements the test of proportional odds. These functions are not documented in detail here as they are primarily for internal use within \pkg{repolr}. There use outside of this setting is not recommended. 
}

\value{
The function \code{summary.repolr} is used to obtain and print a summary of the fitted model.

The fitted model is an object of class \dQuote{\code{repolr}} and has the following values:
                      
  \item{poly.mod}{polynomial model for cut points: a \code{list} with elements
           \code{poly}, \code{polycuts} and \code{space}.}
  \item{y}{the response variable.}
  \item{linear.predictors}{a vector of linear predictors.}
  \item{fitted.values}{a vector of the fitted values.}
  \item{coefficients}{a named vector of regression coefficients.}
  \item{robust.var}{the robust (sandwich) variance matrix.}
  \item{naive.var}{the naive variance matrix.}
  \item{alpha}{an estimate of the correlation parameter.}
  \item{convergence}{a logical variable to reported if convergence was achieved.}
  \item{iter}{the number of iterations.}
  \item{grad1}{first derivative of generalized variance at convergence.}
  \item{grad2}{second derivative of generalized variance at convergence.}
  \item{crit}{convergence criterion.}
  \item{po.test}{results of po test: a \code{list} with elements
  \code{po.stat}, \code{po.df} and \code{po.chi}.}
}

\references{
Parsons NR, Costa ML, Achten J, Stallard N. Repeated measures proportional odds logistic regression analysis of ordinal score data in the statistical software package R. \emph{Computational Statistics and Data Analysis} 2009; 53:632-641.\cr

Parsons NR, Edmondson RN, Gilmour SG. A generalized estimating equation method for fitting autocorrelated ordinal score data with an application in horticultural research. \emph{Journal of the Royal Statistical Society C} 2006; 55:507-524.\cr

Stiger TR, Barnhart HX, Williamson JM. Testing proportionality in the proportional odds model fitted with GEE. \emph{Statistics in Medicine} 1999; 18:1419-1433.\cr

Parsons NR. Proportional-odds models for repeated composite and long ordinal outcome scales. \emph{Statistics in Medicine} 2013; 32:3181-3191.\cr

}
\keyword{models}
\seealso{
\code{\link{QIC}}, \code{\link{polycuts}}, \code{\link{work.corr}}
}


\examples{

data(HHSpain)
mod.0 <- repolr(HHSpain~Sex*Time, data=HHSpain, categories=4, subjects="Patient",
            times=c(1,2,5), corr.mod="uniform", alpha=0.5)
summary(mod.0)
summary(update(mod.0, diffmeth = "numeric"))
summary(update(mod.0, fixed = TRUE, alpha =0.5))

}
