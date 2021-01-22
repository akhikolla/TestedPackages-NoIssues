#' @title Blinded Sample Size Reestimation for Longitudinal Count Data using the NB-INAR(1) Model
#' @description \code{bssr.nb.inar1} fits blinded observations and recalculates the sample size required for proving a desired alternative when testing for
#' a rate ratio between two groups unequal to one. See 'Details' for more information.
#'
#' @param alpha level (type I error) to which the hypothesis is tested.
#' @param power  power (1 - type II error) to which an alternative should be proven.
#' @param delta the rate ratio which is to be proven.
#' @param x a matrix or data frame containing count data which is to be fitted. Columns correspond to time points, rows to observations.
#' @param n a vector giving the sample size within the control group and the treatment group, respecitvely.
#' @param k planned sample size allocation factor between groups: see 'Details'.
#'
#' @details
#' When testing for differences between rates \eqn{\mu_C} and \eqn{\mu_T} of two groups, a control and a treatment group respectively, we usually
#' test for the ratio between the two rates, i.e. \eqn{\mu_T/\mu_C = 1}. The ratio of the two rates is refered to as \eqn{\delta}, i.e.
#' \eqn{\delta = \mu_T/\mu_C}.
#'
#' \code{bssr.nb.inar1} gives back the required sample size for the control and treatment group required to prove an existing
#' alternative \code{theta} with a specified power \code{power} when testing the null hypothesis \eqn{H_0: \mu_T/\mu_C \ge 1} to level \code{alpha}.
#' Nuisance parameters are estimated through the blinded observations \code{x}, thus not further required.
#'
#' for sample sizes \eqn{n_C} and \eqn{n_T} of the control and treatment group, respectively, the argument \code{k} is the desired
#' sample size allocation factor at the end of the study, i.e. \eqn{k = n_T/n_C}.
#'
#' @return \code{rnbinom.inar1} returns the required sample size within the control group and treatment group.
#'
#' @source \code{rnbinom.inar1} uses code contributed by Thomas Asendorf.
#'
#' @seealso \code{\link{rnbinom.inar1}} for information on the NB-INAR(1) model, \code{\link{n.nb.inar1}} for calculating
#' initial sample size required when performing inference, \code{\link{fit.nb.inar1}} for calculating
#' initial parameters required when performing sample size estimation
#'
#' @examples
#' #Calculate required sample size to find significant difference with
#' #80% probability when testing the Nullhypothesis H_0: mu_T/mu_C >= 1
#' #assuming the true effect delta is 0.8 and rate, size and correlation
#' #parameter in the control group are 2, 1 and 0.5, respectively.
#'
#' estimate<-n.nb.inar1(alpha=0.025, power=0.8, delta=0.8, muC=2, size=1, rho=0.5, tp=7, k=1)
#'
#' #Simulate data
#' placebo<-rnbinom.inar1(n=50, size=1, mu=2, rho=0.5, tp=7)
#' treatment<-rnbinom.inar1(n=50, size=1, mu=1.6, rho=0.5, tp=7)
#'
#' #Blinded sample size reestimation
#' blinded.data<-rbind(placebo, treatment)[sample(1:100),]
#' estimate<-bssr.nb.inar1(alpha=0.025, power=0.8, delta=0.8, x=blinded.data, n=c(50,50), k=1)
#' summary(estimate)
#' @export
bssr.nb.inar1<-function(alpha,power,delta,x, n, k){
  tp<-ncol(x)
  datenNA <- tp-rowSums(is.na(x))
  pars<-optim(c(mean(x, na.rm=T), var(x[,1], use="pairwise.complete.obs")/mean(x[,1], na.rm=T)^2, 0.4), minFuncBlinded, lower=c(10,10,10)^-5, upper=c(10,10,1-10^-5), method="L-BFGS-B",
              daten=x, dataNA = datenNA, n=n, delta=delta)$par

  ssest<-n.nb.inar1(alpha=alpha,power=power,delta=delta,muC=(1+k)/(delta*k+1)*pars[1],size=pars[2],rho=pars[3],tp=tp,k=k)
  result<-list(n=ssest$n, alpha=ssest$alpha, power=ssest$power, delta=ssest$delta, mu=ssest$muC, size=ssest$size, rho=ssest$rho, tp=ssest$tp, k=ssest$k, model=ssest$model)
  class(result)<-"bssrest"
  result
}

