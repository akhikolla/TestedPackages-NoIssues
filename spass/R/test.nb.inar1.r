#' @title Testing Hypotheses in NB-INAR(1) model
#' @description \code{test.nb.inar1} tests hypotheses for rate ratios of two groups in an NB-INAR(1) model
#'
#' @param dataC a matrix or data frame containing count data from the control group. Columns correspond to time points, rows to observations.
#' @param dataE a matrix or data frame containing count data from the experiment group. Columns correspond to time points, rows to observations.
#' @param h0     the value against which h is tested, see 'Details'.
#'
#' @details the function \code{test.nb.inar1} tests for the null hypothesis \eqn{\lambda_T/\lambda_C = h0} against the alternative \eqn{\lambda_T/\lambda_C \neq h_0}.
#' For attaining estimates, method of moments estimators are used.
#'
#' @return \code{test.nb.inar1} returns effect size, standard error, Z-statistic and p-value attained through standard normal approximation.
#'
#' @source \code{test.nb.inar1} uses code contributed by Thomas Asendorf.
#'
#' @seealso \code{\link{rnbinom.inar1}} for information on the NB-INAR(1) model, \code{\link{n.nb.inar1}} for calculating
#' initial sample size required when performing inference, \code{\link{fit.nb.inar1}} for calculating
#' initial parameters required when performing sample size estimation
#'
#' @examples
#' set.seed(8)
#' groupE<-rnbinom.inar1(n=1000, size=0.6, mu=2, rho=0.8, tp=6)
#' groupC<-rnbinom.inar1(n=1000, size=0.6, mu=2, rho=0.8, tp=6)
#'
#' test.nb.inar1(dataC=groupC, dataE=groupE, h0=1)
#'
#' @export

test.nb.inar1<-function(dataC, dataE, h0 = 1){
  xE<-dataE
  xC<-dataC
  hyp<-h0
  n.E<-nrow(xE)
  n.C<-nrow(xC)
  tp<-ncol(xE)

  rho.hat <- 1/(n.E+n.C-2)*((n.E-1)*sum(cor(xE, use="pairwise.complete.obs"), na.rm=T)+(n.C-1)*sum(cor(xC, use="pairwise.complete.obs"), na.rm=T))

  lambda.E <- mean(xE, na.rm=T)
  lambda.C <- mean(xC, na.rm=T)

  varterm.E <- 1/tp*1/(n.E-1)*sum((xE-lambda.E)^2, na.rm=T)/lambda.E^2
  varterm.C <- 1/tp*1/(n.C-1)*sum((xC-lambda.C)^2, na.rm=T)/lambda.C^2

  effect<-(log(lambda.E/lambda.C)-log(hyp))
  stderr<- sqrt(1/tp^2*rho.hat*(varterm.E/n.E+varterm.C/n.C))
  Z<-effect/stderr
  erg<-c(effect, stderr, Z, min(pnorm(Z), 1-pnorm(Z))*2)
  names(erg)<-c("effect", "stderr", "Z", "p-value")
  erg
}

