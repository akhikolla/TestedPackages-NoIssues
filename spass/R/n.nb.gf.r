#' @title Sample Size Calculation for Comparing Two Groups when observing Longitudinal Count Data with marginal Negative Binomial Distribution and underlying Gamma Frailty with Autoregressive Correlation Structure of Order One
#' @description \code{n.nb.gf} calculates required sample sizes for testing trend parameters in a Gamma frailty model
#'
#'
#' @param alpha  level (type I error) to which the hypothesis is tested.
#' @param power  power (1 - type II error) to which an alternative should be proven.
#' @param lambda the set of trend parameters assumed to be true at the beginning prior to trial onset
#' @param size   dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer (see \code{\link{rnbinom.gf}}).
#' @param rho    correlation coefficient of the autoregressive correlation structure of the underlying Gamma frailty. Must be between 0 and 1 (see \code{\link{rnbinom.gf}}).
#' @param tp     number of observed time points. (see \code{\link{rnbinom.gf}})
#' @param k      sample size allocation factor between groups: see 'Details'.
#' @param h      hypothesis to be tested. The function must return a single value when evaluated on lambda.
#' @param hgrad  gradient of function h
#' @param h0     the value against which h is tested, see 'Details'.
#' @param trend  the trend which assumed to underlying in the data.
#' @param approx numer of iterations in numerical calculation of the sandwich estimator, see 'Details'.
#'
#' @details
#' The function calculates required samples sizes for testing trend parameters of trends in longitudinal negative binomial data. The underlying
#' one-sided null-hypothesis is defined by \eqn{H_0: h(\eta, \lambda) \geq h_0} vs. the alternative \eqn{H_A: h(\eta, \lambda) < h_0}. For testing
#' these hypothesis, the program therefore requires a function \code{h} and a value \code{h0}.
#'
#' \code{n.nb.gf} gives back the required sample size for the control and treatment group, to prove an existing alternative \eqn{h(\eta, \lambda) - h_0}
#' with a power of \code{power} when testing at level \code{alpha}. For sample sizes \eqn{n_C} and \eqn{n_T} of the control and treatment group, respectively, the argument \code{k} is the
#' sample size allocation factor, i.e. \eqn{k = n_T/n_C}.
#'
#' When calculating the expected sandwich estimator required for the sample size, certain terms can not be computed analytically and have
#' to be approximated numerically. The value \code{approx} defines how close the approximation is to the true expected sandwich estimator.
#' High values of \code{approx} provide better approximations but are compuationally more expensive.
#'
#' @return \code{n.nb.gf} returns the required sample size within the control group and treatment group.
#'
#' @source \code{n.nb.gf} uses code contributed by Thomas Asendorf.
#'
#' @seealso \code{\link{rnbinom.gf}} for information on the Gamma frailty model, \code{\link{fit.nb.gf}} for calculating
#' initial parameters required when performing sample size estimation, \code{\link{bssr.nb.gf}} for blinded
#' sample size reestimation within a running trial.
#'
#'
#' @examples
#' ##The example is commented as it may take longer than 10 seconds to run.
#' ##Please uncomment prior to execution.
#'
#' ##Example for constant rates
#' #h<-function(lambda.eta){
#' #   lambda.eta[2]
#' #}
#' #hgrad<-function(lambda.eta){
#' #   c(0, 1, 0)
#' #}
#'
#' ##We assume the rate in the control group to be exp(lambda[1]) = exp(0) and an
#' ##effect of lambda[2] = -0.3. The \code{size} is assumed to be 1 and the correlation
#' ##coefficient \code{\rho} 0.5. At the end of the study, we would like to test
#' ##the treatment effect specified in lambda[2], and therefore define function
#' ##\code{h} and value \code{h0} accordingly.
#'
#' #estimate<-n.nb.gf(lambda=c(0,-0.3), size=1, rho=1, tp=6, k=1, h=h, hgrad=hgrad,
#' #   h0=0.2, trend="constant", approx=20)
#' #summary(estimate)
#'
#' ##Example for exponential trend
#' #h<-function(lambda.eta){
#' #   lambda.eta[3]
#' #}
#' #hgrad<-function(lambda.eta){
#' #   c(0, 0, 1, 0)
#' #}
#'
#' #estimate<-n.nb.gf(lambda=c(0, 0, -0.3/6), size=1, rho=0.5, tp=7, k=1, h=h, hgrad=hgrad,
#' #   h0=0, trend="exponential", approx=20)
#' #summary(estimate)
#'
#' @export

n.nb.gf <- function(alpha=0.025, power=0.8, lambda, size, rho, tp, k=1, h, hgrad, h0, trend=c("constant", "exponential", "custom"), approx=20){

  type <- which(trend[1] == c("constant", "exponential", "custom"))

  if(type == 3){
    stop("Custom trends are not yet supported in this version.")
  }

  J<-mlFirstJExp(c(lambda,size), rho, k, tp, type, approx)
  H<-mlFirstHExp(c(lambda,size), k, tp, type)

  sigma.2<-(hgrad(c(lambda, size))%*%solve(H)%*%J%*%solve(H)%*%hgrad(c(lambda,size)))[1]

  nC<-(qnorm(1-power)+qnorm(alpha))^2*sigma.2/((1+k)*(h(c(lambda, size))-h0)^2)
  nE<-k*nC
  erg<-c(nC, nE)
  names(erg)<-c("Control", "Treatment")

  result<-list(n = erg, alpha=alpha, power=power, lambda=lambda, size=size, rho=rho, tp=tp, k=k, h=h, hgrad=hgrad, h0=h0, trend=trend, approx=approx, model="GF")
  class(result)<-"ssest"
  return(result)
}
