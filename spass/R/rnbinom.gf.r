#' @title Generate Time Series with Negative Binomial Distribution and Multivariate Gamma Frailty with Autoregressive Correlation Structure of Order One
#' @description \code{rnbinom.gf} generates one or more independent time series following the Gamma frailty model. The generated data has negative binomial marginal distribution and the underlying multivariate Gamma frailty an autoregressive covariance structure.
#'
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param size  dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#' @param mu    vector of means of time points: see 'Details'.
#' @param rho   correlation coefficient of the underlying autoregressive Gamma frailty. Must be between 0 and 1.
#' @param tp    number of observed time points.
#'
#' @details
#' The generated marginal negative binomial distribution with mean \code{mu} = \eqn{\mu} and \code{size} = \eqn{\eta} has density
#' \deqn{(\mu/(\mu+\eta))^x \Gamma(x + \eta)/(\Gamma(x+1)\Gamma(\eta)) (\eta/(\mu+\eta))^\eta}
#' for \eqn{0 < \mu}, \eqn{0 < \eta} and \eqn{x=0, 1, 2, ...}. Hereby, each entry of vector \code{mu} corresponds to
#' one time point. Therefore, each timepoint can have its distinct mean.
#'
#' Within the Gamma frailty model, the correlation between two frailties of time points \eqn{t} and \eqn{s} for \code{rho} = \eqn{\rho} is given by
#' \deqn{\rho^|t-s|}
#' for \eqn{0 \le \rho \le 1}. Note: this does not correspond to the correlation of observations.
#'
#' @return \code{rnbinom.gf} returns a matrix of dimension \code{n} x \code{tp} with marginal negative binomial
#' distribution with means \code{mu}, common dispersion parameter \code{size} and a correlation induce by the autoregressive
#' multivariate Gamma frailty.
#'
#' @source \code{rnbinom.gf} computes observations from a Gamma frailty model by \emph{Fiocco et. al. 2009} using code contributed by Thomas Asendorf.
#'
#' @references Fiocco M, Putter H, Van Houwelingen JC, (2009), A new serially correlated gamma-frailty process for longitudinal count data \emph{Biostatistics} Vol. 10, No. 2, pp. 245-257.
#'
#' @examples
#' set.seed(8)
#' random<-rnbinom.gf(n=1000, size=0.6, mu=1:6, rho=0.8, tp=6)
#' cor(random)
#'
#' #Check the marginal distribution of time point 3
#' plot(table(random[,3])/1000, xlab="Probability", ylab="Observation")
#' lines(0:26, dnbinom(0:26, mu=3, size=0.6), col="red")
#' legend("topright",legend=c("Theoretical Marginal Distribution", "Observed Distribution"),
#'   col=c("red", "black"), lty=1, lwd=c(1,2))
#'
#' @export
rnbinom.gf <- function(n, size, mu, rho, tp){
  Z <- matrix(0, nrow=n, ncol=tp)

  if(length(n)>1){
    n<-length(n)
  }
  if(length(tp)>1){
    tp<-length(tp)
  }
  n<-round(n)
  tp<-round(tp)
  if(size <= 0){
    stop("dispersion parameter must be greater than 0")
  }
  if(any(mu < 0)){
    stop("mean parameter must be greater or equal to 0")
  }
  if(rho < 0 | rho > 1){
    stop("correlation coefficient must be between 0 and 1")
  }
  if(tp <= 0){
    stop("number of time points must be at least 1")
  }
  if(n <= 0){
    stop("number of observations must be at least 1")
  }

  power.ij<-matrix(0, nrow=tp, ncol=tp)
  for(i in 1:tp){
    for(j in 1:tp){
      power.ij[i,j]<-abs(i-j)
    }
  }

  if(n != 1){
    X.i.plus<-matrix(rgamma(tp*n, shape=size*(1-rho)*rho^c(tp:1), rate=size), ncol=tp, byrow=T)
    X.plus.j<-matrix(rgamma(tp*n, shape=size*(1-rho)*rho^c(1:tp), rate=size), ncol=tp, byrow=T)
    X.plus.plus<-matrix(rgamma(n, shape=size*rho^(tp+1), rate=size), ncol=1)
    X.i.j<-array(rgamma(n*tp^2, shape=size*(1-rho)^2*rho^c(power.ij), size), dim=c(tp, tp, n))

    for(t in 1:tp){
      if(t==1 | t==tp){
        Z[,t]<-rowSums(as.matrix(X.i.plus[,1:t]))+rowSums(as.matrix(X.plus.j[,t:tp]))+X.plus.plus+colSums(X.i.j[1:t,t:tp,])
      }else{
        Z[,t]<-rowSums(as.matrix(X.i.plus[,1:t]))+rowSums(as.matrix(X.plus.j[,t:tp]))+X.plus.plus+apply(X.i.j[1:t,t:tp,],3,sum)
      }
    }
  }else{
    X.i.plus<-matrix(rgamma(tp, shape=size*(1-rho)*rho^c(tp:1), rate=size), ncol=tp, byrow=T)
    X.plus.j<-matrix(rgamma(tp, shape=size*(1-rho)*rho^c(1:tp), rate=size), ncol=tp, byrow=T)
    X.plus.plus<-matrix(rgamma(1, shape=size*rho^(tp+1), rate=size), ncol=1)
    X.i.j<-matrix(rgamma(tp^2, shape=size*(1-rho)^2*rho^c(power.ij), size), ncol=tp, nrow=tp)

    for(t in 1:tp){
      Z[,t]<-sum(as.matrix(X.i.plus[,1:t])) + sum(as.matrix(X.plus.j[,t:tp])) + X.plus.plus + sum(X.i.j[1:t, t:tp])
    }
  }

  mu.matrix <- matrix(mu, ncol = tp, nrow = n, byrow = T)
  Y<-matrix(rpois(n*tp, lambda=Z*mu.matrix), nrow=n)
  Y
}


