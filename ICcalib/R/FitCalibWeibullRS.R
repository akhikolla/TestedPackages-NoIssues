#### FitCalibWeibullRS function
## The function takes the following
## w - a matrix. Each row is observation and each column is questionnaire time in the interval. w equal to Inf once
# an observation is censore/had the event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questionnaire result for all participents.
###
# The function returns the estimates of Weibull scale and shape paramters for interval-censored time to event data.

#' @title Fitting Weibull Risk-Set Calibration Models
#' @description Fits Weibull risk-set calibration models for time-to-exposure from interval-censored data. The exposure is a binary covariate measured
#' in intermittent times. This function fits a calibration model at each main event time point, using only members of the risk set at that time point.
#' @param w A matrix of time points when measurements on the binary covariate were obtained.
#' @param w.res A matrix of measurement results of the binary covariate. Each measurement corresponds to the time points in \code{w}
#' @param tm Vector of observed main event time or censoring time.
#' @param event Vector of censoring indicators. \code{1} for event \code{0} for censored
#' @param lower  A value to replace zero in the left point of the interval, Default: 1e-04
#' @param upper A value to replace infinity in the right point of the interval, Default: 200
#' @return A 2-column matrix with the shape and scale parameter for each time-point at which a calibration model was fitted. 
#' @details In case of an error in the model-fitting at a certain time point, a Weibull calibration 
#' model is fitted and used for that time point.
#' @examples 
#' # Simulate data set
#' sim.data <- ICcalib:::SimCoxIntervalCensSingle(n.sample = 200, lambda = 0.1, 
#'                                                alpha = 0.25, beta0 = log(0.5), 
#'                                                mu = 0.2, n.points = 2, 
#'                                                weib.shape = 1, weib.scale = 2)
#' # Fit Weibull risk-set calibration models for the conditional covariate 
#' # starting-time distributions
#' ICcalib::FitCalibWeibullRS(w = sim.data$w, w.res = sim.data$w.res, 
#'                            tm = sim.data$obs.tm, event = sim.data$delta)
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @seealso 
#'  \code{\link[fitdistrplus]{fitdistcens}}, \code{\link[ICcalib]{FitCalibWeibull}}
#' @rdname FitCalibWeibullRS
#' @export 
#' @importFrom fitdistrplus fitdistcens
FitCalibWeibullRS <- function(w, w.res, tm, event, lower = 0.0001, upper = 200)
{
r <- sum(event)
n.fail <- 0
event.index <- which(event==1)
lr.for.fit.all <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
weib.params <- matrix(nrow = r, ncol = 2)
for (j in 1:r)
{
  point <- tm[event.index[j]]
  lr.for.fit <- lr.for.fit.all[tm>=point, ]  # Keep only observations in the risk set
  lr.for.fit <- lr.for.fit[!(lr.for.fit[, 1]==0 & lr.for.fit[, 2]==Inf),]
  colnames(lr.for.fit) <- c("left","right")
  lr.for.fit[lr.for.fit==Inf] <- upper
  lr.for.fit[lr.for.fit==0] <- lower
  fit.weib <- tryCatch(fitdistrplus::fitdistcens(censdata = lr.for.fit, distr = "weibull"),   error=function(e) {e})
  if (inherits(fit.weib, "error")) { 
    weib.params[j,] <- FitCalibWeibull(w, w.res)
    warning(paste("In point", point, "Calibration was used instead of risk set calibration"))
    n.fail <- n.fail + 1
    } else if (fit.weib$estimate[1] > 20 | fit.weib$estimate[2] < 1/1000)
      {
      weib.params[j,] <- FitCalibWeibull(w, w.res)
      warning(paste("In point", point, "Calibration was used instead of risk set calibration"))
      n.fail <- n.fail + 1
      } else   {    
        weib.params[j,] <- fit.weib$estimate
      }
}
if (n.fail > 0) {warning(paste("In ", round(100*n.fail/r,0), "% of the event times there were no sufficient data to fit a risk-set calibration model"))}
if (n.fail/r > 0.5) stop("In more of 50% of the event times there were no sufficient data to fit a risk-set calibration model")
return(weib.params)
}

