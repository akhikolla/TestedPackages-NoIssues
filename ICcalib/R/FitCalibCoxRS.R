#' @title Fitting Proportional Hazards Risk-Set Calibration Models with Covariates
#' @description \code{FitCalibCoxRS} fits proportional hazards risk-set calibration models for time-to-exposure from interval-censored data with covariates. The exposure is a binary covariate measured
#' in intermittent times. The covariates (\code{Q}) are associated with the time-to-exposure. This function fits a calibration model at each main event time point,
#' using only members of the risk set at that time point.
# @param w A matrix of time points when measurements on the binary covariate were obtained.
# @param w.res A matrix of measurement results of the binary covariate. Each measurement corresponds to the time points in \code{w}
# @param Q Matrix of covariates for PH calibration model
# @param hz.times Times used for calculating the baseline hazard function of a PH calibration model
# @param tm Vector of observed main event time or censoring time
# @param n.int The number of interior knots to be used, see \code{ICsurv::fast.PH.ICsurv.EM}, Default: 5
# @param order the order of the basis functions. See \code{ICsurv::fast.PH.ICsurv.EM}, Default: 2
# @param event Vector of censoring indicators. \code{1} for event \code{0} for censored
# @return A list of Cox PH model fits, each supplemented with the knots and order used for the I-splines.
# @details In case of an error in the model-fitting at a certain time point, a proportional hazards calibration 
#' model is fitted (for all the data) and used for that time point.
# @examples 
# sim.data <- ICcalib:::SimCoxIntervalCensCox(n.sample = 50, lambda = 0.1, alpha = 0.25, 
# beta0 = log(0.2), gamma.q = c(log(0.75), log(2.5)), 
# gamma.z = log(1.5), mu = 0.2, n.points = 2)
# # The baseline hazard for the calibration model is calculated in observation times
# cox.hz.times <- sort(unique(sim.data$tm)) 
# # Fit proprtional hazards calibration model
# FitCalibCoxRS(w = sim.data$w, w.res = sim.data$w.res, Q = sim.data$Q, hz.times = cox.hz.times, 
#               tm = sim.data$tm, event = sim.data$delta, n.int = 5, order = 1)
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @seealso 
#'  \code{\link[ICsurv]{fast.PH.ICsurv.EM}},  \code{\link[ICcalib]{FitCalibCox}}
#' @rdname FitCalibCoxRSInts
#' @export 
#' @importFrom ICsurv fast.PH.ICsurv.EM
FitCalibCoxRS <- function(w, w.res, Q, hz.times, tm, n.int = 5, order = 2 , event)
{
r <- sum(event)
n.fail <- 0
event.index <- which(event==1)
lr.for.fit.all <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
Q.all <- Q
all.fit.cox.res <- list()
for (j in 1:r)
{
  point <- tm[event.index[j]]
  # Keep only observations in the risk set
  lr.for.fit <- lr.for.fit.all[tm>=point,]  
  Q <- Q.all[tm>=point,]
  # Take out noninformative observations
  Q <- Q[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),]
  lr.for.fit <- lr.for.fit[!(lr.for.fit[,1]==0 & lr.for.fit[,2]==Inf),]
  #
  colnames(lr.for.fit) <- c("left","right")
  d1 <- lr.for.fit[,1]==0
  d3 <- lr.for.fit[,2]==Inf
  d2 <- 1 - d1 - d3
  fit.cox.point <- tryCatch(ICsurv::fast.PH.ICsurv.EM(d1 = d1, d2 = d2, d3 = d3,Li = lr.for.fit[,1],
                                        Ri = lr.for.fit[,2], n.int = n.int, order = order,  Xp = Q, g0 =rep(1,n.int + order), b0 = rep(0,ncol(Q)),
                                        t.seq = hz.times, tol = 0.001), error = function(e){e})
  while(inherits(fit.cox.point, "error") & n.int >= 2) { 
    n.int <- n.int - 1
    fit.cox.point <- tryCatch(ICsurv::fast.PH.ICsurv.EM(d1 = d1, d2 = d2, d3 = d3,Li = lr.for.fit[,1],
                                          Ri = lr.for.fit[,2], n.int = n.int, order = order,  Xp = Q, g0 =rep(1,n.int + order), b0 = rep(0,ncol(Q)),
                                          t.seq = hz.times, tol = 0.001), error = function(e){e})
  }
  if (n.int<2) {
    fit.cox.point <- FitCalibCox(w = w, w.res = w.res, Q = Q, hz.times = hz.times, n.int = n.int, order = order)
    warning(paste("In point", point, "Calibration was used instead of risk set calibration"))
    n.fail <- n.fail + 1
  }   else {
    ti <- c(lr.for.fit[d1 == 0,1], lr.for.fit[d3 == 0,2])
    fit.cox.point$knots <-   seq(min(ti) - 1e-05,  max(ti) + 1e-05, length.out = (n.int + 2))
    fit.cox.point$order <- order
  }
  all.fit.cox.res[[j]] <- fit.cox.point
}
if (n.fail > 0) {warning(paste("In ", round(100*n.fail/r,0), "% of the event times there were no sufficient data to fit a risk-set calibration model"))}
if (n.fail/r > 0.5) stop("In more of 50% of the event times there were no sufficient data to fit a risk-set calibration model")
return(all.fit.cox.res)
}
