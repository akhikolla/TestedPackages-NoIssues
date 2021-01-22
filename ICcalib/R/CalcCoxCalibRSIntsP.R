#### CalcCoxCalibRSIntsP function
## The function takes the following
## w - a matrix. Each row is observation and each column is questionnaire time in the interval. w equal to Inf once
# an observation is censore/had the event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questionnaire result for all participents.
## point - scalar. The time of the risk set in the main analysis. In terms of the paper, t.
## fit.cox.rs.ints - the result of FitCalibCoxRSInts. This used for the actual risk-set calibration
###
# The function returns a vector with individual predictions for P(X(t)=1|history(time t)). 
# For observations with X(a(t))=1 the above probability is 1 by definition and this is what the
# function returns for them.
#### The following package is needed: fitdistrplus
#### The following functions are used: CalcAuxatPoint (R function)
#' @title Calculating the probabilities of positive binary exposure status at a given time point using proportional hazards grouped risk-set
#'  calibration models
#' @description For a given time point, calculate the probability of positive exposure value  for multiple observations (participants). 
#' The function uses the results of proportional hazards grouped risk-set calibration model fit, and given covariates and collected data on the history 
#' of the binary exposure for each participant. 
#' @param w A matrix of time points when measurements on the binary covariate were obtained.
#' @param w.res A matrix of measurement results of the binary covariate. Each measurement corresponds to the time points in \code{w}
#' @param point The time point at which the probabilities are estimated
#' @param fit.cox.rs.ints The result of \code{FitCalibCoxRSInts} on the interval-censored data
#' @param hz.times Times used for calculating the baseline hazard function from PH calibration model
#' @param Q Matrix of covariates for the PH calibration model
#' @param pts.for.ints Points defining the intervals for grouping risk-sets (first one has to be zero). Should be sorted from zero up
#' @return A vector of estimated probabilities of positive exposure status at time \code{point}.
# @details DETAILS
#' @examples 
#' set.seed(17)
#' sim.data <- ICcalib:::SimCoxIntervalCensCox(n.sample = 100, lambda = 0.1, 
#'                                             alpha = 0.25, beta0 = 0, 
#'                                             gamma.q = c(log(0.75), log(2.5)), 
#'                                             gamma.z = log(1.5), mu = 0.2, 
#'                                             n.points = 2)
#' # The baseline hazard for the calibration model is calculated in observation times
#' cox.hz.times <- sort(unique(sim.data$obs.tm)) 
#' # Fit proprtional hazards calibration model
#' fit.cox.rs.ints <- FitCalibCoxRSInts(w = sim.data$w, w.res = sim.data$w.res, 
#'                                      Q = sim.data$Q, hz.times = cox.hz.times, 
#'                                      n.int = 5, order = 2, pts.for.ints = seq(0,4,1), 
#'                                      tm = sim.data$obs.tm, event = sim.data$delta)
#' # Calculate the conditional probabilities of binary covariate=1 at time one
#' probs <- CalcCoxCalibRSIntsP(w = sim.data$w, w.res = sim.data$w.res, point = 1,
#'                              fit.cox.rs.ints = fit.cox.rs.ints,
#'                              pts.for.ints = seq(0,4,1), Q = sim.data$Q, 
#'                              hz.times = cox.hz.times)
#' summary(probs)
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname CalcCoxCalibRSIntsP
#' @export 
CalcCoxCalibRSIntsP <- function(w, w.res, point, fit.cox.rs.ints, hz.times,  Q, pts.for.ints)
{
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  fit.cox.int <- fit.cox.rs.ints[[findInterval(point, pts.for.ints)]]
  hz <- fit.cox.int$hz
  Qb <- Q%*%fit.cox.int$b
  ## Calculate hazard for the point, first baseline hazard, then add covariates:
  interval.point <- FindIntervalCPP(point = point, w = t(as.matrix(hz.times)))
  if (interval.point==1) {base.hz.point <- hz[1]*point/hz.times[1] } else {
    if(interval.point==length(hz.times)+1) {base.hz.point <- hz[length(hz.times)]
    } else {
      base.hz.point <- hz[interval.point-1] +  (hz[interval.point]-hz[interval.point-1])*(point-hz.times[interval.point-1])/
                    (hz.times[interval.point]-hz.times[interval.point-1])#Extrapolation
    }}
  prob.at.point <- 1-exp(-base.hz.point*exp(Qb[p.point==0,]))
  prob.at.a.point <- 1-CalcSurvFromCox(fit.cox = fit.cox.int,Qb = Qb[p.point==0,], points = a.point[p.point==0], hz.times = hz.times)
  p.point[p.point==0] <- (prob.at.point - prob.at.a.point)/(1 - prob.at.a.point)
  return(p.point)
}
