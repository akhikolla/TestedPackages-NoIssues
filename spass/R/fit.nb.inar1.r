#' @title Fitting Longitudinal Data with Negative Binomial Marginal Distribution and Autoregressive Correlation Structure of Order One: NB-INAR(1)
#' @description \code{fit.nb.inar1} fits data using the maximum likelihood of a reparametrized NB-INAR(1) model.
#'
#' @param x a matrix or data frame containing count data which is to be fitted. Columns correspond to time points, rows to observations.
#' @param lower vector of lower bounds for estimated parameters \code{mu}, \code{size} and \code{rho}, respectively.
#' @param upper vector of upper bounds for estimated parameters \code{mu}, \code{size} and \code{rho}, respectively.
#' @param method algorithm used for minimization of the likelihood, see \code{\link{optim}} for details.
#' @param start vector of starting values for estimated parameters \code{mu}, \code{size} and \code{rho}, respectively, used for optimization.
#'
#' @details the function \code{fit.nb.inar1} fits a reparametrization of the NB-INAR(1) model as found in McKenzie (1986). The reparametrized model
#' assumes equal means and dispersion parameter between time points with an autoregressive correlation structure. The function is especially useful
#' for estimating parameters for an initial sample size calculation using \code{\link{n.nb.inar1}}. The fitting function allows for incomplete follow up,
#' but not for intermittent missingness.
#'
#'
#' @return \code{fit.nb.inar1} return estimates of the mean \code{mu}, dispersion parameter \code{size} and correlation coefficient \code{rho}.
#'
#' @source \code{fit.nb.inar1} uses code contributed by Thomas Asendorf.
#'
#' @seealso \code{\link{rnbinom.inar1}} for information on the NB-INAR(1) model, \code{\link{n.nb.inar1}} for calculating
#' initial sample size required when performing inference, \code{\link{bssr.nb.inar1}} for blinded
#' sample size reestimation within a running trial, \code{\link{optim}} for more information on the used minimization algorithms.
#'
#' @references McKenzie Ed (1986), Autoregressive Moving-Average Processes with Negative-Binomial and Geometric Marginal Distributions. \emph{Advances in Applied Probability} Vol. 18, No. 3, pp. 679-705.
#'
#' @examples
#' #Generate data from the NB-INAR(1) model
#' set.seed(8)
#' random<-rnbinom.inar1(n=1000, size=1.5, mu=2, rho=0.6, tp=7)
#'
#' estimate<-fit.nb.inar1(random)
#' estimate
#' @export

# just testing!

fit.nb.inar1<-function(x, lower=rep(10, 3)^-5, upper=c(10^5, 10^5, 1-10^-5), method="L-BFGS-B", start){
  if(is.vector(x)){
    x<-t(as.matrix(x))
  }else{
    x<-as.matrix(x)
  }
  n<-nrow(x)
  tp<-ncol(x)

  if(n < 10 & tp < 5){
    warning("x contains few observations or few time points. Minimization may be unstable.")
  }
  if(sum(is.na(x)) > 1){
    warning("x contains missing values. fit.nb.inar1 only handles incomplete observations, NOT intermittent missingness.")
  }

  if(is.null(n) | n==1){
    x<-c(x)
    if(missing(start)){
      start<-c(mean(x, na.rm=T), var(x, na.rm=T), max(cor(x[-1], x[-length(x)]), 0.05))
    }
    result<-optim(par=start, fn=minFunc, method=method, lower=lower, upper=upper, daten=x, dataNA=sum(!is.na(x)))$par
  }else{
    if(missing(start)){
      start<-c(mean(x, na.rm=T), mean(diag(var(x, na.rm=T))), max(mean(diag(var(x, na.rm=T)[-1,])), 0.05))
    }
    result<-optim(par=start, fn=minFuncMult, method=method, lower=lower, upper=upper, daten=x, dataNA=rowSums(!is.na(x)), n=n)$par
  }
  names(result)<-c("mu", "size", "rho")
  result
}

