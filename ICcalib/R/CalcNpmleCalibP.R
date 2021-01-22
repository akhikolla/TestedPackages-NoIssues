#### CalcNpmleP function
## The function takes the following
## w - a matrix. Each row is observation and each column is questionnaire time in the interval. w equal to Inf once
# an observation is censore/had the event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questionnaire result for all participents.
## point - scalar. The time of the risk set in the main analysis. In terms of the paper, t.
## fit.npmle - the result of icenReg::ic_np on the whole data. This is used for the actual calibration
###
# The function returns a vector with individual predictions for P(X(t)=1|history(time t)). 
# For observations with X(a(t))=1 the above probability is 1 by definition and this is what the
# function returns for them.
#' @title Calculating the probabilities of positive binary exposure status at a given time point using a nonparametric calibration model 
#' @description For a given time point, calculate the probability of positive exposure value  for multiple observations (participants). 
#' The function uses the results of a nonparametric calibration model fit, and given collected data on the history 
#' of the binary exposure for each participant. 
#' @param w A matrix of time points when measurements on the binary covariate were obtained.
#' @param w.res A matrix of measurement results of the binary covariate. Each measurement corresponds to the time points in \code{w}
#' @param point The time point at which the probabilities are estimated.
#' @param fit.npmle The result of \code{icenReg::ic_np} on the interval-censored data
#' @return A vector of estimated probabilities of positive exposure status at time \code{point}.
# @details DETAILS
#' @examples 
#' sim.data <- ICcalib:::SimCoxIntervalCensSingle(n.sample = 200, lambda = 0.1, 
#'                                                alpha = 0.25, beta0 = log(0.5), 
#'                                                mu = 0.2, n.points = 2, 
#'                                                weib.shape = 1, weib.scale = 2)
#' # Fit nonparametric calibration model
#' fit.npmle <- FitCalibNpmle(w = sim.data$w, w.res = sim.data$w.res)
#' # Calculate the conditional probabilities of binary covariate=1 at time one
#' probs <- CalcNpmleCalibP(w = sim.data$w, w.res = sim.data$w.res, 
#'                          point = 1, fit.npmle = fit.npmle)
#' summary(probs)
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname CalcNpmleCalibP
#' @export 
CalcNpmleCalibP <- function(w, w.res, point, fit.npmle)
{
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  prob.at.point <- 1-CalcSurvFromNPMLE(probs = fit.npmle$p_hat, Tbull = fit.npmle$T_bull_Intervals,
                                       points = point)
  prob.at.a.point <- 1-CalcSurvFromNPMLE(probs = fit.npmle$p_hat, Tbull = fit.npmle$T_bull_Intervals,
                                         points = a.point[p.point==0])
  p.point[p.point==0] <- (prob.at.point - prob.at.a.point)/(1-prob.at.a.point)
  return(p.point)
}
