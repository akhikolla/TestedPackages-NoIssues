#' @title Sample Size Calculation for Comparing Two Groups when observing Longitudinal Count Data with marginal Negative Binomial Distribution and Autoregressive Correlation Structure of Order One: NB-INAR(1)
#' @description \code{n.nb.inar1} calculates the required sample size for proving a desired alternative when testing for
#' a rate ratio between two groups unequal to one. Also gives back power for a specified sample size. See 'Details' for more information.
#'
#' @param alpha level (type I error) to which the hypothesis is tested.
#' @param power  power (1 - type II error) to which an alternative should be proven.
#' @param delta the rate ratio which is to be proven.
#' @param muC   the rate observed within the control group.
#' @param size  dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer (see \code{\link{rnbinom.inar1}}).
#' @param rho   correlation coefficient of the underlying autoregressive correlation structure. Must be between 0 and 1 (see \code{\link{rnbinom.inar1}}).
#' @param tp    number of observed time points. (see \code{\link{rnbinom.inar1}})
#' @param k     sample size allocation factor between groups: see 'Details'.
#' @param npow  sample size for which a power is to be calculated. Can not be specified if power is also specified.
#' @param nmax  maximum total sample size of both groups. If maximum is reached a warning message is broadcasted.
#'
#' @details
#' When testing for differences between rates \eqn{\mu_C} and \eqn{\mu_T} of two groups, a control and a treatment group respectively, we usually
#' test for the ratio between the two rates, i.e. \eqn{\mu_T/\mu_C = 1}. The ratio of the two rates is refered to as \eqn{\delta}, i.e.
#' \eqn{\delta = \mu_T/\mu_C}.
#'
#' \code{n.nb.inar1} gives back the required sample size for the control and treatment group required to prove an existing
#' alternative \code{theta} with a specified power \code{power} when testing the null hypothesis \eqn{H_0: \mu_T/\mu_C \ge 1} to level \code{alpha}.
#' If \code{power} is not specified but instead \code{npow}, the power achieved with a total sample size of \code{npow} is calculated.
#'
#' For sample sizes \eqn{n_C} and \eqn{n_T} of the control and treatment group, respectively, the argument \code{k} is the
#' sample size allocation factor, i.e. \eqn{k = n_T/n_C}.
#'
#'
#' @return \code{rnbinom.inar1} returns the required sample size within the control group and treatment group.
#'
#' @source \code{rnbinom.inar1} uses code contributed by Thomas Asendorf.
#'
#' @seealso \code{\link{rnbinom.inar1}} for information on the NB-INAR(1) model, \code{\link{fit.nb.inar1}} for calculating
#' initial parameters required when performing sample size estimation, \code{\link{bssr.nb.inar1}} for blinded
#' sample size reestimation within a running trial.
#'
#' @examples
#' #Calculate required sample size to find significant difference with
#' #80% probability when testing the Nullhypothesis H_0: mu_T/mu_C >= 1
#' #assuming the true effect delta is 0.8 and rate, size and correlation
#' #parameter in the control group are 2, 1 and 0.5, respectively.
#'
#' estimate<-n.nb.inar1(alpha=0.025, power=0.8, delta=0.8, muC=2, size=1, rho=0.5, tp=7, k=1)
#' summary(estimate)
#'
#' estimate<-n.nb.inar1(alpha=0.025, npow=200, delta=0.8, muC=2, size=1, rho=0.5, tp=7, k=1)
#' summary(estimate)
#' @export

n.nb.inar1<-function(alpha,power=NULL,delta,muC,size,rho,tp,k,npow=NULL,nmax=Inf){
  beta<-power
  if((is.null(power) & is.null(npow)) | (!is.null(power) & !is.null(npow))){
    stop("either power or npow have to be specified")
  }
  if(rho == 1){
    dep.par<-tp^2
  }else{
    dep.par<-2/(rho-1)^2*(rho^(tp+1)-(rho-1)*(tp+1)-1)-tp
  }
  lambdaquer<-1/(1+k)*muC*(delta*k+1)

  if(is.null(npow)){
    n<-(qnorm(alpha)+qnorm(1-beta))^2*dep.par/(tp^2*log(delta)^2)*((1+k*delta)^2/((1+k)*k*delta*lambdaquer)+1/size*(1+1/k))
    ns<-round(c(n, k*n))
    pow<-1-pnorm(-sqrt(ns[1]/((1+k*delta)^2/((1+k)*k*delta*lambdaquer)+1/size*(1+1/k))*(tp^2*log(delta)^2)/dep.par)-qnorm(alpha))
  }
  if(is.null(power)){
    pow<-1-pnorm(-sqrt(round(npow/(1+k))/((1+k*delta)^2/((1+k)*k*delta*lambdaquer)+1/size*(1+1/k))*(tp^2*log(delta)^2)/dep.par)-qnorm(alpha))
    ns<-c(round(npow/(1+k)), npow-round(npow/(1+k)))
  }

  if(sum(ns)>=nmax){
    pow.nmax<-1-pnorm(-sqrt(round(nmax/(1+k))/((1+k*delta)^2/((1+k)*k*delta*lambdaquer)+1/size*(1+1/k))*(tp^2*log(delta)^2)/dep.par)-qnorm(alpha))
    warning(paste("calculated sample size exceeds nmax =", nmax, "\n expected power using nmax subjects:", round(pow.nmax,2)))
  }

  names(ns)<-c("Control", "Treatment")
  result<-list(n = ns, alpha=alpha, power=pow, delta=delta, muC=muC, size=size, rho=rho, tp=tp, k=k, model="NB-INAR(1)")
  class(result)<-"ssest"
  return(result)
}

