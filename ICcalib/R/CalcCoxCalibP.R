#' @title Calculating the probabilities of positive binary exposure status at a given time point using a proportional hazards calibration model 
#' @description For a given time point, calculate the probability of positive exposure value  for multiple observations (participants). 
#' The function uses the results of a proportional hazards calibration model fit, and given covariates and collected data on the history 
#' of the binary exposure for each participant. 
#' @param w A matrix of time points when measurements on the binary covariate were obtained.
#' @param w.res A matrix of measurement results of the binary covariate. Each measurement corresponds to the time points in \code{w}
#' @param point The time point at which the probabilities are estimated
#' @param fit.cox The result of \code{icenReg::ic_sp} on the interval-censored data
#' @param hz.times Times used for calculating the baseline hazard function from PH calibration model
#' @param Q Matrix of covariates for the PH calibration model
#' @return A vector of estimated probabilities of positive exposure status at time \code{point}.
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
#' fit.cox <- FitCalibCox(w = sim.data$w, w.res = sim.data$w.res, Q = sim.data$Q, 
#'                        hz.times = cox.hz.times, n.int = 5, order = 2)
#' # Calculate the conditional probabilities of binary covariate=1 at time one
#' probs <- CalcCoxCalibP(w = sim.data$w, w.res = sim.data$w.res, point = 1,
#'                        Q = sim.data$Q, fit.cox = fit.cox, hz.times = cox.hz.times)
#' summary(probs)
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname CalcCoxCalibP
#' @export 
CalcCoxCalibP <- function(w, w.res, point, fit.cox, hz.times,  Q)
{
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  hz <- fit.cox$hz
  Qb <- Q%*%fit.cox$b
  ## Calculate hazard for the point, first baseline hazard, then add covariates:
  interval.point <- FindIntervalCPP(point = point, w = t(as.matrix(hz.times)))
  if (interval.point==1) {base.hz.point <- hz[1]*point/hz.times[1] } else {
    if(interval.point==length(hz.times)+1) {base.hz.point <- hz[length(hz.times)]
    } else {
      base.hz.point <- hz[interval.point-1] +  (hz[interval.point]-hz[interval.point-1])*(point-hz.times[interval.point-1])/
                    (hz.times[interval.point]-hz.times[interval.point-1]) #Extrapolation
    }}
  prob.at.point <- 1-exp(-base.hz.point*exp(Qb[p.point==0,]))
  prob.at.a.point <- 1-CalcSurvFromCox(fit.cox = fit.cox,Qb = Qb[p.point==0,], points = a.point[p.point==0], hz.times = hz.times)
  p.point[p.point==0] <- (prob.at.point - prob.at.a.point)/(1-prob.at.a.point)
  return(p.point)
}
