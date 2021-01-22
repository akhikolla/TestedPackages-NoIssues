## point - scalar. The time of the risk set in the main analysis. In terms of the paper, t.
## weib.params - the shape and scale parameters from the Weibull calibration fitting to the interval-cenosed time to exposure/treatment
###
# The function returns a vector with individual predictions for P(X(t)=1|history(time t)). 
# For observations with X(a(t))=1 the above probability is 1 by definition and this is what the
# function returns for them.
# The function calculates prediction for all observations, even though predictions for observations outside 
# the risk set are not used
#### The following functions is used: CalcAuxatPoint (R function)
#' @title Calculating the probabilities of positive binary exposure status at a given time point using a Weibull calibration model 
#' @description For a given time point, calculate the probability of positive exposure value  for multiple observations (participants). 
#' The function uses the results of a Weibull calibration model fit, and given collected data on the history 
#' of the binary exposure for each participant. 
#' @param w A matrix of time points when measurements on the binary covariate were obtained.
#' @param w.res A matrix of measurement results of the binary covariate. Each measurement corresponds to the time points in \code{w}
#' @param point The time point at which the probabilities are estimated.
#' @param weib.params A bivariate vector. Shape and scale parameters of the Weibull calibration model.  
#' @return A vector of estimated probabilities of positive exposure status at time \code{point}.
# @details DETAILS
#' @examples 
#' # Simulate data set
#' sim.data <- ICcalib:::SimCoxIntervalCensSingle(n.sample = 200, lambda = 0.1, 
#'                                                alpha = 0.25, beta0 = log(0.5), 
#'                                                mu = 0.2, n.points = 2, 
#'                                                weib.shape = 1, weib.scale = 2)
#' # Fit a Weibull calibration model for the covariate starting time distribution
#' calib.weib.params <- FitCalibWeibull(w = sim.data$w, w.res = sim.data$w.res)
#' # Calculate the conditional probabilities of binary covariate=1 at time one
#' probs <- CalcWeibullCalibP(w = sim.data$w, w.res = sim.data$w.res, point = 1,
#'                            weib.params = calib.weib.params)
#' summary(probs)
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @seealso 
#'  \code{\link[stats]{Weibull}}
#' @rdname CalcWeibullCalibP
#' @export 
#' @importFrom stats pweibull
CalcWeibullCalibP <- function(w, w.res, point, weib.params)
{
  weib.shape <- weib.params[1]
  weib.scale <- weib.params[2]
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  prob.at.point <- stats::pweibull(point, shape = weib.shape,scale = weib.scale)
  prob.at.a.point <- stats::pweibull(a.point[p.point==0], shape = weib.shape, scale = weib.scale)
  p.point[p.point==0] <- (prob.at.point - prob.at.a.point)/(1-prob.at.a.point)
  return(p.point)
}
