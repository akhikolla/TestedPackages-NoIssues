#### CalcNpmleRSP function
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
#' @title Calculating the probabilities of positive binary exposure status at a given time point using a nonparametric risk-set calibration models 
#' @description For a given time point, calculate the probability of positive exposure value  for multiple observations (participants). 
#' The function first fits the nonparametric risk-set calibration models at each main event time point and then calculates the probabilities
#' of positive binary exposure status.
#' @param w A matrix of time points when measurements on the binary covariate were obtained.
#' @param w.res A matrix of measurement results of the binary covariate. Each measurement corresponds to the time points in \code{w}
#' @param point The time point at which the probabilities are estimated
#' @param obs.tm Vector of observed main event time or censoring time
#' @return A vector of estimated probabilities of positive exposure status at time \code{point}.
#' @details This function calculates the NPMLE at each main event time point and then provides the estimated probabilities for positive 
#' exposure status at time \code{point}.
#' @examples 
#' # Simulate data set
#' sim.data <- ICcalib:::SimCoxIntervalCensSingle(n.sample = 200, lambda = 0.1, 
#'                                                alpha = 0.25, beta0 = log(0.5), 
#'                                                mu = 0.2, n.points = 2, 
#'                                                weib.shape = 1, weib.scale = 2)
#' # Calculate the conditional probabilities of binary covariate=1 at time one
#' # Unlike CalcNpmle, CalcNpmleRSP includes the calibration model fitting
#' probs <- CalcNpmleRSP(w = sim.data$w, w.res = sim.data$w.res, point = 1, 
#'                       obs.tm = sim.data$obs.tm)
#' summary(probs)
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @seealso 
#'  \code{\link[icenReg]{ic_np}}
#' @rdname CalcNpmleRSP
#' @export 
#' @importFrom icenReg ic_np
CalcNpmleRSP <- function(w, w.res, point, obs.tm)
{
  lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
  colnames(lr.for.fit) <- c("left","right")
  lr.for.fit <- lr.for.fit[obs.tm >= point,] # Keep only observations in the risk set
  fit.npmple.rs <-  tryCatch(icenReg::ic_np(cbind(left,right) ~ 0, data = lr.for.fit),   error=function(e) {e})
  if (inherits(fit.npmple.rs, "error")) { 
    #fit.npmle <- FitCalibNpmle(w = w, w.res = w.res)
    #p.point <- CalcNpmleCalibP(w = w, w.res = w.res, point = point, fit.npmle = fit.npmle)
    #warning(paste("In point", point, "Calibration was used instead of risk set calibration"))
    #return(p.point)
    stop(paste("In point", point, "there were no sufficient data to fit a risk-set calibration model"))
    }
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  prob.at.point <- 1-CalcSurvFromNPMLE(probs = fit.npmple.rs$p_hat, Tbull = fit.npmple.rs$T_bull_Intervals,
                                       points = point)
  prob.at.a.point <- 1-CalcSurvFromNPMLE(probs = fit.npmple.rs$p_hat, Tbull = fit.npmple.rs$T_bull_Intervals,
                                         points = a.point[p.point==0])
  # If the support of the estimated distriubtion ends before someones last available questionnire, the we get a contradiction
  # because \hat{F}(a)=1 but X(a)=0. In this rare case, we just do carry forward the zero.
  prob.at.a.point[prob.at.a.point==1] <-  -prob.at.point
  # If the last questionnire result was X(a)=1 then no problem caused above since p.point==1
  p.point[p.point==0] <- (prob.at.point - prob.at.a.point)/(1-prob.at.a.point)
  return(p.point)
}
