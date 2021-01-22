#' @title Testing Hypotheses in Gamma Frailty models
#' @description \code{test.nb.gf} tests hypotheses for certain trends in Gamma frailty models
#'
#' @param dataC a matrix or data frame containing count data from the control group. Columns correspond to time points, rows to observations.
#' @param dataE a matrix or data frame containing count data from the experiment group. Columns correspond to time points, rows to observations.
#' @param h      hypothesis to be tested. The function must return a single value when evaluated on lambda.
#' @param hgrad  gradient of function h
#' @param h0     the value against which h is tested, see 'Details'.
#' @param trend  the trend which assumed to be underlying in the data.
#' @param H0 indicates if the sandwich estimator is calculated under the null hypothesis or alternative.
#' @param one.sided indicates if the hypothesis should be tested one- or two-sided
#' @param ... Arguments to be passed to function \code{fit.nb.gf()}.
#'
#' @details the function \code{test.nb.gf} tests for the null hypothesis \eqn{h(\eta, \lambda) = h_0} against the alternative \eqn{h(\eta, \lambda) \neq h_0}.
#' The fitting function allows for incomplete follow up, but not for intermittent missingness.
#'
#' If parameter H0 is set to TRUE, the hessian and outer gradient are calculated under the assumption that \code{lambda[2]} \eqn{\geq} \code{h0} if
#' \code{trend = "constant"} or \code{lambda[3]} \eqn{\geq} \code{h0} if \code{trend = "exponential"}.
#'
#' @return \code{test.nb.gf} returns effect size, standard error, Z-statistic and p-value attained through standard normal approximation.
#'
#' @source \code{test.nb.gf} uses code contributed by Thomas Asendorf.
#'
#' @seealso \code{\link{rnbinom.gf}} for information on the Gamma Frailty model, \code{\link{n.nb.gf}} for calculating
#' initial sample size required when performing inference, \code{\link{fit.nb.gf}} for calculating
#' initial parameters required when performing sample size estimation.
#'
#' @references Fiocco M, Putter H, Van Houwelingen JC, (2009), A new serially correlated gamma-frailty process for longitudinal count data \emph{Biostatistics} Vol. 10, No. 2, pp. 245-257.
#'
#' @examples
#' #Create data from two groups
#' random<-get.groups(n=c(100,100), size=c(0.7, 0.7), lambda=c(0.8, 0), rho=c(0.6, 0.6),
#'   tp=7, trend="constant")
#'
#' #Define hypothesis
#' h<-function(lambda.eta){
#'   lambda.eta[2]
#' }
#' hgrad<-function(lambda.eta){
#'   c(0, 1, 0)
#' }
#' test.nb.gf(dataC=random[101:200,], dataE=random[1:100,], h=h, hgrad=hgrad, h0=0,
#'   trend="constant", H0=FALSE)
#' @export

test.nb.gf <- function(dataC, dataE, h, hgrad, h0=0, trend=c("constant", "exponential", "custom"), H0=FALSE, one.sided=TRUE, ...){

  if(trend == "custom"){
    stop("Custom trends are not yet supported in this version.")
  }

  groupC <- dataC
  groupE <- dataE

  n<-c(nrow(groupC), nrow(groupE))

  erg.fit<-fit.nb.gf(dataC=groupC, dataE=groupE, trend=trend, method="L-BFGS-B", H0=H0, h0=h0, ...)

  H<-erg.fit$hessian
  lambda<-erg.fit$par
  J<-erg.fit$ogradient

  effect <- (h(lambda)-h0)
  stderr <- sqrt((hgrad(lambda)%*%(solve(H)%*%J%*%solve(H))%*%hgrad(lambda))[1,1])/sqrt(sum(n))
  Z<-effect/stderr

  if(one.sided == FALSE){
    erg<-c(effect, stderr, Z, min(pnorm(Z), 1-pnorm(Z))*2)
  }else{
    erg<-c(effect, stderr, Z, pnorm(Z))
  }

  names(erg)<-c("effect", "stderr", "Z", "p-value")
  erg
}
