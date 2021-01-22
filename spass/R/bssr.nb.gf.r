#' @title Blinded Sample Size Reestimation for Longitudinal Count Data with marginal Negative Binomial Distribution and underlying Gamma Frailty with Autoregressive Correlation Structure of Order One
#' @description \code{bssr.nb.gf} fits blinded observations and recalculates the sample size required for sustaining power at desired alternative when testing for
#' trend parameters in a Gamma frailty models. See 'Details' for more information.
#'
#' @param data a matrix or data frame containing count data which is to be fitted. Columns correspond to time points, rows to observations.
#' @param alpha level (type I error) to which the hypothesis is tested.
#' @param power  power (1 - type II error) to which an alternative should be proven.
#' @param delta  the relevant effect size, which is assumed to be true, see 'Details'.
#' @param h0     the value against which h is tested, see 'Details'.
#' @param tp     number of observed time points. (see \code{\link{rnbinom.gf}})
#' @param k      sample size allocation factor between groups: see 'Details'.
#' @param trend  the trend which assumed to underlying in the data.
#' @param approx numer of iterations in numerical calculation of the sandwich estimator, see 'Details'.
#'
#' @details
#' The function recalculates a sample size for testing in constant and exponential trends.
#'
#' Under a constant trend, the means in control and experiment group are equal to \eqn{\lambda_1} and  \eqn{\lambda_1 + \lambda_2}, respectively.
#' The treatment effect \code{delta} is therefore equal to \eqn{\lambda_2}.
#'
#' Under an exponential trend, the means in control and experiment group are equal to \eqn{exp(\lambda_1+t \cdot \lambda_2)} and  \eqn{\lambda_1 + t\cdot \lambda_2 + t\cdot \lambda_3}, respectively.
#' The treatment effect \code{delta} is therefore equal to \eqn{\lambda_3}.
#'
#' \code{bssr.nb.gf} returns the required sample size for the control and treatment group required to prove an existing
#' alternative \code{delta} with a specified power \code{power} when testing the null hypothesis \eqn{H_0: \delta \ge h_0} at level \code{alpha}.
#' Nuisance parameters are estimated through the blinded observations \code{data}, thus not further required.
#' For sample sizes \eqn{n_C} and \eqn{n_T} of the control and treatment group, respectively, the argument \code{k} is the desired
#' sample size allocation factor at the end of the study, i.e. \eqn{k = n_T/n_C}.
#'
#' @return \code{bssr.nb.gf} returns the required sample size within the control group and treatment group.
#'
#' @source \code{bssr.nb.gf} uses code contributed by Thomas Asendorf.
#'
#' @seealso \code{\link{rnbinom.gf}} for information on the Gamma Frailty model, \code{\link{n.nb.gf}} for calculating
#' initial sample size required when performing inference, \code{\link{fit.nb.gf}} for calculating
#' initial parameters required when performing sample size estimation.
#'
#' @examples
#' ##The example is commented as it may take longer than 10 seconds to run.
#' ##Please uncomment prior to execution.
#'
#' ##Example for constant rates
#' #set.seed(12)
#' #h<-function(lambda.eta){
#' #   lambda.eta[2]
#' #}
#' #hgrad<-function(lambda.eta){
#' #   c(0, 1, 0)
#' #}
#'
#' ##Calculate initial sample size
#' #estimate<-n.nb.gf(lambda=c(0,-0.3), size=1, rho=0.5, tp=6, k=1, h=h, hgrad=hgrad,
#' #   h0=0, trend="constant", approx=20)
#'
#' ##Generate and permutate data with different nuisance parameters
#' #random<-get.groups(n=round(estimate$n/2), size=c(0.8, 0.8), lambda=c(0.5, -0.3),
#' #   rho=c(0.4, 0.4), tp=6, trend="constant")
#' #random<-random[sample(1:nrow(random), nrow(random)), ]
#'
#' ##Recalculate sample size with data
#' #reestimate<-bssr.nb.gf(data=random, alpha=0.025, power=0.8, delta=-0.3, h0=0,
#' #   tp=6, k=1, trend="constant", approx = 20)
#'
#' #summary(reestimate)
#'
#' @import stats
#' @export
bssr.nb.gf<-function(data, alpha=0.025, power=0.8, delta, h0=0, tp, k, trend=c("constant", "exponential", "custom"), approx = 20){

  data<-as.matrix(data)

  if(trend == "constant"){
    type <- 1
    start <- c(log(mean(data, na.rm=T)), mean(data, na.rm=T)^2/(var(c(data), na.rm=T)-mean(data, na.rm=T)))
    lower <- c(-Inf, 1e-10);upper <- c(Inf,Inf)
    h<-function(lambda.eta){
      lambda.eta[2]
    }
    hgrad<-function(lambda.eta){
      c(0, 1, 0)
    }
  }else if(trend == "exponential"){
    type <- 2
    start <- c(log(mean(data[,1], na.rm=T)), log(mean(data[,2], na.rm=T))-log(mean(data[,1], na.rm=T)), max(0.01, mean(data, na.rm=T)^2/(var(c(data), na.rm=T)-mean(data, na.rm=T))))
    lower <- c(-Inf, -Inf, 1e-10); upper <- c(Inf, Inf, Inf)
    h<-function(lambda.eta){
      lambda.eta[3]
    }
    hgrad<-function(lambda.eta){
      c(0, 0, 1, 0)
    }

  }else if(trend == "custom"){
    stop("Custom trends are not yet supported in this version.")
  }

  erg <- optim(start, fn=mlFirstBlinded, lower = lower, upper=upper, method = "L-BFGS-B", group=data, n=nrow(data), tp=rowSums(!is.na(data)),
               type=type, theta=delta, k=k)

  if(trend == "constant"){
    y <- c(erg$par[1], delta, erg$par[length(erg$par)])
  }else if(trend == "exponential"){
    y <- c(erg$par[1:2], delta, erg$par[length(erg$par)])
  }

  rho <- optim(cor(data[,1], data[,2], use="pairwise.complete.obs"), y=y, fn=mlSecondBlinded, lower = 0, upper=1, method = "Brent", group=data, n=nrow(data),
               type=type, tp = rowSums(!is.na(data)), k=k)$par

  result<-n.nb.gf(lambda=y[1:(length(y)-1)], size=y[length(y)], rho=rho, h=h, h0=h0, hgrad=hgrad,
          trend=trend, tp=tp, alpha=alpha, power=power, k=k, approx=approx)
  class(result)<-"bssrest"
  result$model <- "GF"
  result$delta <- delta
  result
}


