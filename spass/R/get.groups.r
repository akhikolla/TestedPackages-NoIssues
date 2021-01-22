#' @title Generate Time Series with Negative Binomial Distribution and Multivariate Gamma Frailty with Autoregressive Correlation Structure of Order One with Trend
#' @description \code{rnbinom.gf} generates one or more independent time series following the Gamma frailty model. The generated data has negative binomial marginal distribution and the underlying multivariate Gamma frailty an autoregressive covariance structure.
#'
#' @param n number of observations.
#' @param size  dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#' @param lambda    vector of means of trend parameters.
#' @param rho   correlation coefficient of the underlying autoregressive Gamma frailty. Must be between 0 and 1.
#' @param tp    number of observed time points.
#' @param trend  a string giving the trend which is to be simulated.
#'
#' @details
#' The function relies on \code{\link{rnbinom.gf}} for creating data with underlying constant or exponential trends.
#'
#' @return \code{get.groups} returns a matrix of dimension \code{n} x \code{tp} with marginal negative binomial
#' distribution with means corresponding to trend parameters \code{lambda}, common dispersion parameter \code{size} and a correlation induce by \code{rho},
#' the correlation coefficient of the autoregressive multivariate Gamma frailty.
#'
#' @source \code{rnbinom.gf} computes observations from a Gamma frailty model by \emph{Fiocco et. al. 2009} using code contributed by Thomas Asendorf.
#'
#' @references Fiocco M, Putter H, Van Houwelingen JC, (2009), A new serially correlated gamma-frailty process for longitudinal count data \emph{Biostatistics} Vol. 10, No. 2, pp. 245-257.
#'
#' @seealso \code{\link{rnbinom.gf}} for information on the Gamma frailty model.
#'
#' @examples
#' random<-get.groups(n=c(1000,1000), size=c(0.5, 0.5), lambda=c(1, 2), rho=c(0.6, 0.6), tp=7,
#'   trend="constant")
#' head(random)
#'
#' @export
get.groups <- function(n, size, lambda, rho, tp, trend){
  if(trend == "constant"){
    mu1<-rep(exp(lambda[1]+lambda[2]), tp)
    mu2<-rep(exp(lambda[1]), tp)
  }else if(trend == "exponential"){
    mu1<-exp(lambda[1]+(lambda[2]+lambda[3])*0:(tp-1))
    mu2<-exp(lambda[1]+(lambda[2])*0:(tp-1))
  }
  group.E<-rnbinom.gf(n[1], size[1], mu1, rho[1], tp)
  group.C<-rnbinom.gf(n[2], size[2], mu2, rho[2], tp)
  return(rbind(group.E,group.C))
}


