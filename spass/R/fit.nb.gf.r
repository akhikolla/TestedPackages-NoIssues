#' @title Fitting Longitudinal Data from a Gamma Frailty Model with Frailty of Autoregressive Correlation Structure of Order One
#' @description \code{fit.nb.gf} fits data using the pseudo maximum likelihood of a Gamma frailty model
#'
#' @param dataC a matrix containing count data from the control group, which is to be fitted. Columns correspond to time points, rows to observations.
#' @param dataE a matrix containing count data from the experiment group, which is to be fitted. Columns correspond to time points, rows to observations.
#' @param trend  the trend which assumed to underlying in the data.
#' @param lower vector of lower bounds for estimated parameters \code{lambda}, \code{size} and \code{rho}, respectively.
#' @param upper vector of upper bounds for estimated parameters \code{lambda}, \code{size} and \code{rho}, respectively.
#' @param method algorithm used for minimization of the likelihood, see \code{\link{optim}} for details.
#' @param start vector of starting values for estimated parameters \code{mu}, \code{size} and \code{rho}, respectively, used for optimization.
#' @param approx numer of iterations in numerical calculation of the sandwich estimator, see 'Details'.
#' @param rho indicates whether or not to calculate the correlation coefficient of Gamma frailties. Must be TRUE or FALSE.
#' @param H0 indicates whether or not to calculate the hessian and outer gradient matrix under the null hypothesis, see 'Details'.
#' @param h0 the value against which is tested under the null
#'
#' @details the function \code{fit.nb.gf} fits a Gamma frailty model as found in Fiocco (2009). The fitting function allows for incomplete follow up,
#' but not for intermittent missingness.
#'
#' When calculating the expected sandwich estimator required for the sample size, certain terms can not be computed analytically and have
#' to be approximated numerically. The value \code{approx} defines how close the approximation is to the true expected sandwich estimator.
#' High values of \code{approx} provide better approximations but are compuationally more expensive.
#'
#' If parameter H0 is set to TRUE, the hessian and outer gradient are calculated under the assumption that \code{lambda[2]} \eqn{\geq} \code{h0} if
#' \code{trend = "constant"} or \code{lambda[3]} \eqn{\geq} \code{h0} if \code{trend = "exponential"}.
#'
#' @return \code{fit.nb.gf} returns estimates of the trend parameters \code{lambda}, dispersion parameter \code{size},
#' Hessian matrix \code{hessian}, outer gradient product matrix \code{ogradient} and, if inquired, correlation coefficient \code{rho}.
#'
#' @source \code{fit.nb.gf} uses code contributed by Thomas Asendorf.
#'
#' @seealso \code{\link{rnbinom.gf}} for information on the Gamma frailty model, \code{\link{n.nb.gf}} for calculating
#' initial sample size required when performing inference, \code{\link{bssr.nb.gf}} for blinded
#' sample size reestimation within a running trial, \code{\link{optim}} for more information on the used minimization algorithms.
#'
#' @references Fiocco M, Putter H, Van Houwelingen JC, (2009), A new serially correlated gamma-frailty process for longitudinal count data \emph{Biostatistics} Vol. 10, No. 2, pp. 245-257.
#'
#' @examples
#' #Generate data from the Gamma frailty model
#' random<-get.groups(n=c(1000,1000), size=c(0.7, 0.7), lambda=c(0.8, -0.5), rho=c(0.6, 0.6),
#'   tp=7, trend="constant")
#' fit.nb.gf(dataC=random[1001:2000,], dataE=random[1:1000,], trend="constant")
#' @export

fit.nb.gf <- function(dataC, dataE, trend=c("constant", "exponential"),
                      lower, upper, method="L-BFGS-B", start, approx=20, rho = FALSE, H0 = FALSE, h0=0){
  groupC<-as.matrix(dataC)
  groupE<-as.matrix(dataE)

  nC<-nrow(groupC)
  nE<-nrow(groupE)

  tpC<-rowSums(!is.na(groupC))
  tpE<-rowSums(!is.na(groupE))

  tp<-ncol(groupC)

  if(H0==TRUE){
    if(trend == "constant"){
      if(missing(lower)){
        lower.re <- c(-Inf, h0, 1e-10)
        lower <- c(-Inf, -Inf, 1e-10)
      }
      if(missing(upper)){
        upper <- c(Inf, Inf, Inf)
      }
      if(missing(start)){
        start <- c(log(mean(groupE, na.rm=T)), log(mean(groupC, na.rm=T)), max(0.01, mean(groupC, na.rm=T)^2/(var(groupC[,1], na.rm=T)-mean(groupC, na.rm=T))))
        start.re <- c(log(mean(groupE, na.rm=T)), max(log(mean(groupC, na.rm=T)), h0), max(0.01, mean(groupC, na.rm=T)^2/(var(groupC[,1], na.rm=T)-mean(groupC, na.rm=T))))
      }
      erg <- optim(start, mlFirst, gr=mlFirstGrad, lower = lower, upper=upper, method = "L-BFGS-B",
                   groupE=groupE, groupC=groupC, tpE=tpE, tpC=tpC,
                   nE=nE, nC=nC, type=1)
      erg.re <- optim(start.re, mlFirst, gr=mlFirstGrad, lower = lower.re, upper=upper, method = "L-BFGS-B",
                      groupE=groupE, groupC=groupC, tpE=tpE, tpC=tpC,
                      nE=nE, nC=nC, type=1)
      y <- c(erg.re$par[1], erg.re$par[2], erg.re$par[length(erg.re$par)])
      rho <- optim(cor(groupC[,1], groupC[,2], use="pairwise.complete.obs"), y=y, fn=mlSecond, lower = 0, upper=1, method = "Brent", groupE=groupE, groupC=groupC,
                   nE=nE, nC=nC, tpE=tpE, tpC=tpC, type=1)$par

      erg$rho<-rho
      erg$hessian<-mlFirstHExp(y, nE/nC, tp, type=1)
      erg$ogradient<-mlFirstJExp(y, rho, nE/nC, tp, type=1, approx=20)

      names(erg$par) <- c("lambda1", "lambda2", "size")
    }else if(trend == "exponential"){
      if(missing(lower)){
        lower <- c(-Inf, -Inf, -Inf, 1e-10)
        lower.re <- c(-Inf, -Inf, h0, 1e-10)
      }
      if(missing(upper)){
        upper <- c(Inf, Inf, Inf, Inf)
      }
      if(missing(start)){
        start <- c(log(mean(groupC[,1], na.rm=T)), log(mean(groupE[,2], na.rm=T))-log(mean(groupE[,1], na.rm=T)), log(mean(groupC[,2], na.rm=T))-log(mean(groupC[,1], na.rm=T))-log(mean(groupE[,2], na.rm=T))+log(mean(groupE[,1])), max(0.01, mean(groupC, na.rm=T)^2/(var(groupC[,1], na.rm=T)-mean(groupC, na.rm=T))))
        start.re <- c(log(mean(groupC[,1], na.rm=T)), log(mean(groupE[,2], na.rm=T))-log(mean(groupE[,1], na.rm=T)), max(log(mean(groupC[,2], na.rm=T))-log(mean(groupC[,1], na.rm=T))-log(mean(groupE[,2], na.rm=T))+log(mean(groupE[,1], na.rm=T)),h0), max(0.01, mean(groupC, na.rm=T)^2/(var(groupC[,1], na.rm=T)-mean(groupC, na.rm=T))))
      }
      erg <- optim(start, mlFirst, gr=mlFirstGrad, lower = lower, upper=upper, method = "L-BFGS-B",
                   groupE=groupE, groupC=groupC, tpE=tpE, tpC=tpC,
                   nE=nE, nC=nC, type=2)
      erg.re <- optim(start.re, mlFirst, gr=mlFirstGrad, lower = lower, upper=upper, method = "L-BFGS-B",
                   groupE=groupE, groupC=groupC, tpE=tpE, tpC=tpC,
                   nE=nE, nC=nC, type=2)

      y <- c(erg.re$par[1], erg.re$par[2], erg.re$par[3], erg.re$par[length(erg.re$par)])
      rho <- optim(cor(groupC[,1], groupC[,2], use="pairwise.complete.obs"), y=y, fn=mlSecond, lower = 0, upper=1, method = "Brent", groupE=groupE, groupC=groupC,
                   nE=nE, nC=nC, tpE=tpE, tpC=tpC, type=2)$par
      erg$rho<-rho
      erg$hessian<-mlFirstHExp(y, nE/nC, tp, type=2)
      erg$ogradient<-mlFirstJExp(y, rho, nE/nC, tp, type=2, approx=20)

      names(erg$par) <- c("lambda1", "lambda2", "lambda3", "size")
    }
  }else{
    if(trend == "constant"){
      type <- 1
      if(missing(lower)){
        lower <- c(-Inf, -Inf, 1e-10)
      }
      if(missing(upper)){
        upper <- c(Inf, Inf, Inf)
      }
      if(missing(start)){
        start <- c(log(mean(groupE)), log(mean(groupC)), max(0.01, mean(groupC)^2/(var(groupC[,1])-mean(groupC))))
      }
      erg <- optim(start, mlFirst, gr=mlFirstGrad, lower = lower, upper=upper, method = "L-BFGS-B",
                   groupE=groupE, groupC=groupC, tpE=tpE, tpC=tpC,
                   nE=nE, nC=nC, type=1)

      erg$hessian<-mlFirstHObs(erg$par, groupE, groupC, nE, nC, tpE, tpC, type=1)
      erg$ogradient<-mlFirstJObs(erg$par, groupE, groupC, nE, nC, tpE, tpC, type=1)

      names(erg$par) <- c("lambda1", "lambda2", "size")
    }else if(trend == "exponential"){
      type <- 2
      if(missing(lower)){
        lower <- c(-Inf, -Inf, -Inf, 1e-10)
      }
      if(missing(upper)){
        upper <- c(Inf, Inf, Inf, Inf)
      }
      if(missing(start)){
        start <- c(log(mean(groupC[,1])), log(mean(groupE[,2]))-log(mean(groupE[,1])), log(mean(groupC[,2]))-log(mean(groupC[,1]))-log(mean(groupE[,2]))+log(mean(groupE[,1])), max(0.01, mean(groupC)^2/(var(groupC[,1])-mean(groupC))))
      }
      erg <- optim(start, mlFirst, gr=mlFirstGrad, lower = lower, upper=upper, method = "L-BFGS-B",
                   groupE=groupE, groupC=groupC, tpE=tpE, tpC=tpC,
                   nE=nE, nC=nC, type=2)

      erg$hessian<-mlFirstHObs(erg$par, groupE, groupC, nE, nC, tpE, tpC, type=2)
      erg$ogradient<-mlFirstJObs(erg$par, groupE, groupC, nE, nC, tpE, tpC, type=2)

      names(erg$par) <- c("lambda1", "lambda2", "lambda3", "size")
    }

    if(rho==TRUE){
      rho <- optim(cor(dataC[,1], dataC[,2], use="pairwise.complete.obs"), y=erg$par, groupE = as.matrix(dataE), groupC = as.matrix(dataC), nE = nE, nC = nC,
                   tpE = tpE, tpC = tpC, type = type, fn=mlSecond, lower = 0, upper=1, method = "Brent")$par
      erg$rho <- rho
    }
  }
  erg
}


