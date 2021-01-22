#### CalcWeibullRSP function
## The function takes the following
## w - a matrix. Each row is observation and each column is questionnaire time in the interval. w equal to Inf once
# an observation is censore/had the event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questionnaire result for all participents.
## point - scalar. The time of the risk set in the main analysis. In terms of the paper, t.
## obs.tm - the vector of observed times for all observaions. In terms of the paper, T. This is used for finding the risk set.
###
# The function returns a vector with individual predictions for P(X(t)=1|history(time t)). 
# For observations with X(a(t))=1 the above probability is 1 by definition and this is what the
# function returns for them.
# The function calculate prediction for all observations, even though predictions for observations outside 
# the risk set are not used 
#### The following functions is used: CalcAuxatPoint (R function), FindIntervalCalibCPP (cpp function)
#' @title Calculating the probabilities of positive binary exposure status at a given time point using risk-set Weibull calibration models 
#' @description For a given time point, calculate the probability of positive exposure value  for multiple observations (participants). 
#' The function uses the results of a Weibull calibration model fit, and given collected data on the history 
#' of the binary exposure for each participant. 
#' @param w A matrix of time points when measurements on the binary covariate were obtained.
#' @param w.res A matrix of measurement results of the binary covariate. Each measurement corresponds to the time points in \code{w}
#' @param point The time point at which the probabilities are estimated.
#' @param weib.params A bivariate vector. Shape and scale parameters of the Weibull calibration model.  
#' @return A vector of estimated probabilities of positive exposure status at time \code{point}.
#' @details At its present form this function is identical to \code{CalcWeibullCalibP}. This is because the current version of the \code{ICcalib} package
#' (Version 1.0.005), the user loop over the main event times. Then, at each event time point, the user should include the appropriate Weibull
#' parameters as estimated by \code{FitCalibWeibullRS}. 
#' @examples 
#' # Simulate data set
#' sim.data <- ICcalib:::SimCoxIntervalCensSingle(n.sample = 200, lambda = 0.1, 
#'                                                alpha = 0.25, beta0 = log(0.5), 
#'                                                mu = 0.2, n.points = 2, 
#'                                                weib.shape = 1, weib.scale = 2)
#' case.times <- sim.data$obs.tm[sim.data$delta==1]
#' # Fit Weibull risk-set calibration models
#' calib.weib.params <- FitCalibWeibullRS(w = sim.data$w, w.res = sim.data$w.res, 
#'                                        tm = sim.data$obs.tm, 
#'                                        event = sim.data$delta)
#' # Calculate the conditional probabilities of binary covariate=1 at time one
#' probs <- CalcWeibullRSP(w = sim.data$w, w.res = sim.data$w.res, point = 1,
#'                         weib.params = calib.weib.params)
#' summary(probs)
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[stats]{Weibull}}
#' @rdname CalcWeibullRSP
#' @export 
#' @importFrom stats pweibull
CalcWeibullRSP <- function(w, w.res, point, weib.params)
{
  weib.shape <- weib.params[1]
  weib.scale <- weib.params[2]
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  prob.at.point <- stats::pweibull(point, shape = weib.shape, scale = weib.scale)
  prob.at.a.point <- stats::pweibull(a.point[p.point==0], shape = weib.shape, scale = weib.scale)
  p.point[p.point==0] <- (prob.at.point - prob.at.a.point)/(1-prob.at.a.point)
  return(p.point)
}
