#### Functions for fitting ordinery calibration ###
# On Sep 1, 2016, this file included three calibrations: Weibull, Nonparameteric and Cox. 
# The Weibull calibrations return the shape and scale paramters
# The Cox and NP calibrations returns fitted objects.
# Sep 1, 2016: Only the Cox calibration function allows for covariates
################################################################################################################
################### Weibull ###############################################################################################

#' @title Fitting Weibull Calibration Models
#' @description Fits a Weibull calibration model for time-to-exposure from interval-censored data. The exposure is a binary covariate measured
#' in intermittent times.
#' @param w A matrix of time points when measurements on the binary covariate were obtained.
#' @param w.res A matrix of measurement results of the binary covariate. It corresponds to the time points in \code{w}
#' @return A bivariate vector: the estimated Weibull shape and scale parameters
# @details DETAILS
#' @examples
#' # Simulate data set
#' sim.data <- ICcalib:::SimCoxIntervalCensSingle(n.sample = 200, lambda = 0.1, 
#'                                                alpha = 0.25, beta0 = log(0.5), 
#'                                                mu = 0.2, n.points = 2, 
#'                                                weib.shape = 1, weib.scale = 2)
#' # Fit a Weibull calibration model for the covariate starting time distribution
#' ICcalib::FitCalibWeibull(w = sim.data$w, w.res = sim.data$w.res) 
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @seealso 
#'  \code{\link[fitdistrplus]{fitdistcens}}
#' @rdname FitCalibWeibull
#' @export 
#' @importFrom fitdistrplus fitdistcens
#' @importFrom Rcpp compileAttributes
#' @importFrom survival Surv
#' @useDynLib ICcalib, .registration = TRUE
FitCalibWeibull <- function(w,w.res)
{
  lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
  lr.for.fit <- lr.for.fit[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),]
  colnames(lr.for.fit) <- c("left","right")
  lr.for.fit[lr.for.fit==Inf] <- NA
#  lr.for.fit[lr.for.fit==0] <- 0.0001
  #df.cln <- lr.for.fit[apply(lr.for.fit,1,function(x)  sum(is.na(x)))<2,]
  fit.weib <- tryCatch(fitdistrplus::fitdistcens(censdata = lr.for.fit, distr = "weibull"),   error=function(e) {e})
  if (inherits(fit.weib, "error")) { 
    fit.weib <- tryCatch(fitdistrplus::fitdistcens(censdata = lr.for.fit, distr = "weibull", lower = c(0, 0)),   error=function(e) {e})
    if (inherits(fit.weib, "error")) { 
      return(c(NA,NA)) }
  } 
  if (fit.weib$estimate[1]> 20 | fit.weib$estimate[2] < 1/1000)
  {
    return(c(NA,NA)) 
  }
  return(fit.weib$estimate)
}

################################################################################################################
################### Nonparametric (NPMLE)#####################################################################################
#' @title Fitting Nonparametric Calibration Models
#' @description Fits a Weibull calibration model for time-to-exposure from interval-censored data. The exposure is a binary covariate measured
#' in intermittent times.
#' @param w A matrix of time points when measurements on the binary covariate were obtained.
#' @param w.res A matrix of measurement results of the binary covariate. It corresponds to the time points in \code{w}
#' @return The result of NPMLE model fitting for interval-censored data, as obtained by \code{icenReg::ic_np}
# @details DETAILS
#' @examples 
#' # Simulate data set
#' sim.data <- ICcalib:::SimCoxIntervalCensSingle(n.sample = 200, lambda = 0.1, 
#'                                                alpha = 0.25, beta0 = log(0.5), 
#'                                                mu = 0.2, n.points = 2, 
#'                                                weib.shape = 1, weib.scale = 2)
#' # Calcualte the NPMLE for a nonparametric calibration model for the covariate 
#' # starting-time distribution
#' npmle.fit <- ICcalib::FitCalibNpmle(w = sim.data$w, w.res = sim.data$w.res)
#' # Plot the estimated survival function
#' plot(npmle.fit)
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @seealso 
#'  \code{\link[icenReg]{ic_np}}
#' @rdname FitCalibNpmle
#' @export 
#' @importFrom icenReg ic_np
FitCalibNpmle <- function(w,w.res)
{
  lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
  colnames(lr.for.fit) <- c("left","right")
  fit.npmple <- icenReg::ic_np(cbind(left,right)~0,data = lr.for.fit)
  if (inherits(fit.npmple, "error")) { 
    return(NA) 
  } else {
    return(fit.npmple)
  }
}
################################################################################################################
################### Cox ####################################################################################################
#' @title Fitting Proportional Hazards Calibration Models with Covariates
#' @description Fits a proportional hazards calibration model for time-to-exposure from interval-censored data with covariates. The exposure is a binary covariate measured
#' in intermittent times. The covariates (\code{Q}) are associated with the time-to-exposure.
#' @param w A matrix of time points when measurements on the binary covariate were obtained.
#' @param w.res A matrix of measurement results of the binary covariate. It corresponds to the time points in \code{w}
#' @param Q Matrix of covariates for PH calibration model
#' @param hz.times Times used for calculating the baseline hazard function from PH calibration model
#' @param n.int The number of interior knots to be used, see \code{ICsurv::fast.PH.ICsurv.EM}, Default: 5
#' @param order the order of the basis functions. See \code{ICsurv::fast.PH.ICsurv.EM}, Default: 2
#' @return An object created by ICsurv::fast.PH.ICsurv.EM, with additional variables \code{knots} and \code{order}.
# @details DETAILS
#' @examples 
#' sim.data <- ICcalib:::SimCoxIntervalCensCox(n.sample = 200, lambda = 0.1, 
#'                                             alpha = 0.25, beta0 = 0, 
#'                                             gamma.q = c(log(0.75), log(2.5)), 
#'                                             gamma.z = log(1.5), mu = 0.2, 
#'                                             n.points = 2)
#' # The baseline hazard for the calibration model is calculated in observation times
#' cox.hz.times <- sort(unique(sim.data$obs.tm)) 
#' # Fit proprtional hazards calibration model
#' FitCalibCox(w = sim.data$w, w.res = sim.data$w.res, Q = sim.data$Q, 
#'             hz.times = cox.hz.times, n.int = 5, order = 2)
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @seealso 
#'  \code{\link[ICsurv]{fast.PH.ICsurv.EM}}
#' @rdname FitCalibCox
#' @export 
#' @importFrom ICsurv fast.PH.ICsurv.EM
FitCalibCox <- function(w, w.res, Q, hz.times, n.int = 5, order = 2 )
{
  lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
  Q <- as.matrix(Q[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),])
  lr.for.fit <- lr.for.fit[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),]
  colnames(lr.for.fit) <- c("left","right")
  n.s <- nrow(lr.for.fit)
  d1 <- lr.for.fit[,1]==0
  d3 <- lr.for.fit[,2]==Inf
  d2 <- 1 - d1 - d3
  fit.cox <- tryCatch(ICsurv::fast.PH.ICsurv.EM(d1 = d1, d2 = d2, d3 = d3,Li = lr.for.fit[,1],
                               Ri = lr.for.fit[,2], n.int = n.int, order = order,  Xp = Q, g0 =rep(1,n.int + order), b0 = rep(0,ncol(Q)),
                               t.seq = hz.times, tol = 0.001), error = function(e){e})
  while(inherits(fit.cox, "error") & n.int >= 2) { 
    n.int <- n.int - 1
    fit.cox <- tryCatch(ICsurv::fast.PH.ICsurv.EM(d1 = d1, d2 = d2, d3 = d3,Li = lr.for.fit[,1],
                                 Ri = lr.for.fit[,2], n.int = n.int, order = order,  Xp = Q, g0 =rep(1,n.int + order), b0 = rep(0,ncol(Q)),
                                 t.seq = hz.times, tol = 0.001), error = function(e){e})
  }
  if (n.int<2) {return(NA)}   else {
    ti <- c(lr.for.fit[d1 == 0,1], lr.for.fit[d3 == 0,2])
    fit.cox$knots <-   seq(min(ti) - 1e-05,  max(ti) + 1e-05, length.out = (n.int + 2))
    fit.cox$order <- order
    return(fit.cox)
  }
}
######################################################
