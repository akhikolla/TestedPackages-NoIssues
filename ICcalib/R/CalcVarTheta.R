###
#' @title Variance estimation for the main proportional hazards model
#' @description Estimation of the covariance matrix for the parameters of the main proportional hazards model. 
#' This includes the variance of the binary exposure estimate and the other covariates, if included in the model. 
#' Each function correspond to a different calibration (or risk-set calibration model).
#' @param theta Coefficient vector from main PH model. First coefficient corresponds to X, the rest to Z
#' @param tm Vector of observed main event time or censoring time
#' @param event Vector of censoring indicators. \code{1} for event \code{0} for censored
#' @param Z Additional variables for the main model other than the binary covariate
#' @param Q For PH calibration models: additional covariates
#' @param ps A matrix. Rows are observations, columns are time pointas of the events. The entry at the i-th row and j-column is
#' the conditional probability of positive exposure status for observation i at the j-th event time.
#' @param ps.deriv A matrix. Rows are observations, columns are time points of the events.  The derivative of \code{ps} with 
#' respect to the calibration model parameters
#' @param fit.cox For PH calibration models:  The result of \code{icenReg::ic_sp} on the interval-censored data
#' @return The covariance matrix. The first row and column are for the binary exposure. 
# @details DETAILS
#' @examples
#' # Simulate data set
#' sim.data <- ICcalib:::SimCoxIntervalCensSingle(n.sample = 200, lambda = 0.1, alpha = 0.25, 
#'                                                beta0 = log(0.5), mu = 0.2, n.points = 2, 
#'                                                weib.shape = 1, weib.scale = 2)
#' case.times <- sim.data$obs.tm[sim.data$delta==1]
#' # Fit a Weibull calibration model for the covariate starting time distribution
#' calib.weib.params <- FitCalibWeibull(w = sim.data$w, w.res = sim.data$w.res)
#' px <- t(sapply(case.times, CalcWeibullCalibP, w = sim.data$w, 
#'                w.res =  sim.data$w.res, weib.params = calib.weib.params))
#' # Calculate derivative matrices
#' px.deriv.shape <- t(sapply(case.times, ICcalib:::CalcWeibullCalibPderivShape, 
#' w = sim.data$w, w.res =  sim.data$w.res, weib.params = calib.weib.params))
#' px.deriv.scale <- t(sapply(case.times, ICcalib:::CalcWeibullCalibPderivScale, 
#' w = sim.data$w, w.res = sim.data$ w.res, weib.params = calib.weib.params))
#' # Point estimate 
#' est.weib.calib <- optimize(f = ICcalib:::CoxLogLikX,  tm = sim.data$obs.tm, 
#'                            event = sim.data$delta, ps = px, interval = c(-50,50), 
#'                            maximum = TRUE)$maximum
#' # Variance estimate (no addtional covariates)
#' var.beta.wb <- CalcVarThetaWeib(beta = est.weib.calib, etas = calib.weib.params, 
#'                                 tm = sim.data$obs.tm, event = sim.data$delta, 
#'                                 ps = px, ps.deriv.shape = px.deriv.shape, 
#'                                 ps.deriv.scale =  px.deriv.scale, w = sim.data$w, 
#'                                 w.res = sim.data$w.res)
#'  print(est.weib.calib)
#'  print(var.beta.wb)                               
#  \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
# @seealso 
#  \code{\link[MASS]{ginv}}
#' @rdname CalcVar
#' @export 
#' @importFrom MASS ginv
CalcVarParam <- function(theta,  tm, event, Z, Q, ps, ps.deriv, w, w.res, fit.cox)
{
  n <- length(tm)
  n.theta <- length(theta) 
  eta.b <- fit.cox$b
  eta.g <- fit.cox$g
  hessian.eta <- fit.cox$Hessian
  order <- fit.cox$order
  knots <- fit.cox$knots
  n.eta <- length(eta.b) + length(eta.g)
  nabla.eta.Utheta <- matrix(nrow = n.theta, ncol = n.eta, 0)
  for (i in 1:n.eta)
  {
    nabla.eta.Utheta[1,i] <- CalcNablabeetaUbeta(theta = theta, tm = tm, event = event, ps = ps, Z = Z, psDeriv = t(ps.deriv[,i,]) )
    nabla.eta.Utheta[2:n.theta,i] <- CalcNablabeetaUgamma(theta = theta, tm = tm, event = event, ps = ps, Z = Z, psDeriv = t(ps.deriv[,i,]))
  }
  nabla.eta.Utheta <- nabla.eta.Utheta/n
  ### Prep for calculating gradient of eta
  
  lr.for.fit.raw <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
  Qinter <- as.matrix(Q[!(lr.for.fit.raw[,1]==0 & lr.for.fit.raw[,2]==Inf), ])
  lr.for.fit <- lr.for.fit.raw[!(lr.for.fit.raw[,1]==0 & lr.for.fit.raw[,2]==Inf), ]
  d1 <- lr.for.fit[,1]==0
  d3 <- lr.for.fit[,2]==Inf
  d2 <- 1 - d1 - d3
  
  ######
  grad.eta.pers.mat <- matrix(nrow = n, ncol = n.eta,0)
  grad.eta.pers.mat[!(lr.for.fit.raw[,1]==0 & lr.for.fit.raw[,2]==Inf), ] <- CalcGradEtaPers(d1 = d1, d2 = d2, d3 = d3, Li = lr.for.fit[,1],
                                       Ri = lr.for.fit[,2], knots = knots, order = order, eta.g = eta.g, eta.b = eta.b, Q = Qinter)
  b.mat <- CalcbZ(theta = theta, tm = tm, event = event, ps = ps, Z = Z)
  r.mat <- b.mat - t(nabla.eta.Utheta%*%MASS::ginv(hessian.eta)%*%t(grad.eta.pers.mat))
  MM <- array(dim = c(ncol(r.mat),ncol(r.mat),nrow(r.mat)),0)
  for (j in 1:nrow(r.mat))
  {
    MM[,,j] <- r.mat[j,]%*%t(r.mat[j,])
  } 
  meat <- apply(MM, c(1,2), mean)
  bread <- solve(CoxLogLikHess(theta = theta, tm = tm, event = event, ps = ps, Z = Z))
  v.hat <- n*bread%*%meat%*%bread
  if(max(v.hat)>1) {
    for (j in 1:nrow(b.mat))
    {
      MM[,,j] <- b.mat[j, ]%*%t(b.mat[j, ])
    } 
  }
  meat <- apply(MM,c(1,2),mean)
  v.hat <- n*bread%*%meat%*%bread
  return(v.hat)
}

# @title Variance estimation for main proportional hazards model under proportional hazards grouped risk-set calibration models
# @description Calculates the covariance matrix for the parameters of the main calibration model. This includes the variance. 
# @param theta Coefficient vector from main PH model. First coefficient corresponds to X, the rest to Z.
# @param tm Vector of observed main event time or censoring time.
# @param event Vector of censoring indicators. \code{1} for event \code{0} for censored.
# @param Z Additional variables for the main model other than the binary covariate.
# @param Q Matrix of covariates for PH calibration model.
# @param ps A matrix. Rows are observations, columns are time points of the events. The entry at the i-th row and j-column is
# the conditional probability of positive exposure status for observation i at the j-th event time.
# @param ps.deriv A matrix. Rows are observations, columns are time points of the events.  The derivative of \code{ps} with 
# respect to the calibration model parameters.
#' @param w A matrix of time points when measurements on the binary covariate were obtained.
#' @param w.res A matrix of measurement results of the binary covariate. The measurement corresponding to the time points in \code{w}.
#' @param fit.cox.rs.ints For grouped risk-set PH calibraion: The result of \code{FitCalibCoxRSInts} on the interval-censored data
#' @param pts.for.ints For grouped-risk set PH calibraion: Points defining the intervals for grouping risk-sets (first one has to be zero). Should be sorted from zero up
#' @param n.etas.per.fit For grouped-risk set PH calibraion: A vector. Total number of parameters for each PH calibration fit.
# @return The covariance matrix. The first row and column are for the binary exposure.
# @details DETAILS
# @examples 
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
# @seealso 
#  \code{\link[MASS]{ginv}}
#' @rdname CalcVar
#' @export 
#' @importFrom MASS ginv
CalcVarParamRSInts <- function(theta,  tm, event, Z, Q, ps, ps.deriv, w, w.res,  fit.cox.rs.ints,  pts.for.ints, n.etas.per.fit)
{
  n <- length(tm)
  n.theta <- length(theta) 
  n.fits <- length(fit.cox.rs.ints)
  n.eta <- sum(n.etas.per.fit)
  hessian.eta <- matrix(nrow = n.eta, ncol = n.eta, 0)
  nabla.eta.Utheta <- matrix(nrow = n.theta, ncol = n.eta, 0)
  for (j in 1:n.fits)
  {
    point <- pts.for.ints[j]
    n.pars.ints <- n.etas.per.fit[j]
    fit.temp.ints <- fit.cox.rs.ints[[j]]
    if (j > 1)
    {hessian.eta[(sum(n.etas.per.fit[1:(j-1)]) + 1):(sum(n.etas.per.fit[1:j])), 
                (sum(n.etas.per.fit[1:(j-1)]) + 1):(sum(n.etas.per.fit[1:j]))] <- fit.temp.ints$Hessian
    } else {
      hessian.eta[1:n.etas.per.fit[j], 1:n.etas.per.fit[j]] <- fit.temp.ints$Hessian
    }
    
    in.risk.set <- tm >= point
    tm.ints <- tm[in.risk.set]
    event.ints <- event[in.risk.set]
    ps.ints <- ps[ ,in.risk.set]
    Z.ints <- Z[in.risk.set,]
    ps.deriv.ints <- ps.deriv[in.risk.set,,]
    for (i in 1:n.pars.ints)
    {
      if (j>1) {param.index <- sum(n.etas.per.fit[1:(j-1)]) + i} else {param.index <- i}
      nabla.eta.Utheta[1, param.index] <- CalcNablabeetaUbeta(theta = theta, tm = tm.ints, event = event.ints, ps = ps.ints, Z = Z.ints, 
                                                   psDeriv = t(ps.deriv.ints[, param.index, ]) )
      nabla.eta.Utheta[2:n.theta, param.index] <- CalcNablabeetaUgamma(theta = theta, tm = tm.ints, event = event.ints, ps = ps.ints, 
                                                                       Z = Z.ints, psDeriv = t(ps.deriv.ints[, param.index, ]))  
    }
    
  }
  nabla.eta.Utheta <- nabla.eta.Utheta/n
  ### Prep for calculating gradient of eta
  
  lr.for.fit.raw <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
  Qinter <- as.matrix(Q[!(lr.for.fit.raw[,1]==0 & lr.for.fit.raw[,2]==Inf), ])
  tm.inter <- tm[!(lr.for.fit.raw[,1]==0 & lr.for.fit.raw[,2]==Inf)] 
  lr.for.fit <- lr.for.fit.raw[!(lr.for.fit.raw[,1]==0 & lr.for.fit.raw[,2]==Inf), ]
  d1 <- lr.for.fit[,1]==0
  d3 <- lr.for.fit[,2]==Inf
  d2 <- 1 - d1 - d3
  
  ######
  grad.eta.pers.mat <- matrix(nrow = n, ncol = n.eta,0)
  grad.eta.pers.mat[!(lr.for.fit.raw[,1]==0 & lr.for.fit.raw[,2]==Inf), ] <- CalcGradEtaPersRSInts(d1 = d1, d2 = d2, d3 = d3, 
                                                                                                   Li = lr.for.fit[,1], Ri = lr.for.fit[,2], 
                                                                                                   Q = Qinter, 
                                                                                                   fit.cox.rs.ints = fit.cox.rs.ints,
                                                                                                   pts.for.ints = pts.for.ints, tm = tm.inter, 
                                                                                                   n.etas.per.fit = n.etas.per.fit)
    #CalcGradEtaPersRSInts(d1 = d1, d2 = d2, d3 = d3, Li = lr.for.fit[,1], Ri = lr.for.fit[,2], knots = knots, 
    #order = order, eta.g = eta.g, eta.b = eta.b, Z = Zinter)
  b.mat <- CalcbZ(theta = theta, tm = tm, event = event, ps = ps, Z = Z)
  r.mat <- b.mat - t(nabla.eta.Utheta%*%MASS::ginv(hessian.eta)%*%t(grad.eta.pers.mat))
  MM <- array(dim = c(ncol(r.mat),ncol(r.mat),nrow(r.mat)),0)
  for (j in 1:nrow(r.mat))
  {
    MM[,,j] <- r.mat[j,]%*%t(r.mat[j,])
  } 
  meat <- apply(MM, c(1,2), mean)
  bread <- solve(CoxLogLikHess(theta = theta, tm = tm, event = event, ps = ps, Z = Z))
  v.hat <- n*bread%*%meat%*%bread
  if(max(v.hat)>1) {
    for (j in 1:nrow(b.mat))
    {
      MM[,,j] <- b.mat[j, ]%*%t(b.mat[j, ])
    } 
  }
  meat <- apply(MM,c(1,2),mean)
  v.hat <- n*bread%*%meat%*%bread
  return(v.hat)
}


## Weibull, no extra covariates ###
# @title Variance estimation and confidence intervals for the effect of the exposure under Weibull calibration.
# @description Bootstrap calculations of the variance and confidence interval for the for the log hazard-ratio of the binary exposure
# under a Weibull calibration model.
# @param beta Coefficient of the binary covariate
#' @param etas For Weibull calibration: Shape and scale parameters of the Weibull calibration model.
# @param tm Vector of observed main event time or censoring time
# @param event Vector of censoring indicators. \code{1} for event \code{0} for censored
# @param ps A matrix. Rows are observations, columns are time points of the events. The entry at the i-th row and j-column is
# the conditional probability of positive exposure status for observation i at the j-th event time.
#' @param ps.deriv.shape The derivative of \code{ps} with respect to the shape parameter of the Weibull calibration model.
#' @param ps.deriv.scale The derivative of \code{ps} with respect to the scale parameter of the Weibull calibration model.
# @param w A matrix of time points when measurements on the binary covariate were obtained.
# @param w.res A matrix of measurement results of the binary covariate. Each measurement corresponds to the time points in \code{w}
# @return Variance estimate and possibly confidence interval for the log hazard-ratio of the binary exposure under 
# @details DETAILS
# @examples 
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
# @seealso 
#  \code{\link[numDeriv]{hessian}}
#' @rdname CalcVar
#' @export 
#' @importFrom numDeriv hessian
CalcVarThetaWeib <- function(beta, etas, tm, event, ps, ps.deriv.shape, ps.deriv.scale, w, w.res)
{
  n <- length(tm)
  b.vec <- Calcb(beta = beta, tm = tm, event = event, ps = ps)
  nabla.eta.shape.Ubeta <- CalcUbetabeeta(beta = beta, tm = tm, event = event, ps = ps, psDeriv = ps.deriv.shape)
  nabla.eta.scale.Ubeta <- CalcUbetabeeta(beta = beta, tm = tm, event = event, ps = ps, psDeriv = ps.deriv.scale)
  nabla.eta.Ubeta <- c(nabla.eta.shape.Ubeta, nabla.eta.scale.Ubeta)/n
  hess.etas.l.v <- (numDeriv::hessian(func = ICweibLik, x = etas, w = w, w.res = w.res))
  grad.eta.pers <- ICweibGrad(etas = etas, w = w, w.res = w.res)
  r.vec <- b.vec - nabla.eta.Ubeta%*%solve(hess.etas.l.v)%*%t(grad.eta.pers)
  meat <- mean(r.vec^2)  # since beta is one-dimensional here 
  bread <- myFmyHess(beta, tm, event, ps)/n
  var.beta <- (meat/(bread^2))/n
  return(var.beta)
}
# @title Variance estimation and confidence intervals for the effect of the exposure under Weibull risk-set calibration.
# @description Bootstrap calculations of the variance and confidence interval for the for the log hazard-ratio of the binary exposure
# under  Weibull risk-set calibration models.
#' @param beta Coefficient of the binary covariate. The analogue of theta for non-PH calibration models
#' @param etas.matrix For Weibull risk-set calibration: Two-columns matrix. Each row contains shape and scale parameters from a Weibull risk-set calibration model
# @param tm Vector of observed main event time or censoring time
# @param event Vector of censoring indicators. \code{1} for event \code{0} for censored
#' @param ps.rs A matrix. Rows are observations, columns are time points of the events. The entry at the i-th row and j-column is
#' the conditional probability of positive exposure status for observation i at the j-th event time.
#' @param ps.deriv.shape.rs For Weibull risk-set calibration:The derivative of \code{ps} with respect to the shape parameters of the Weibull risk-set calibration models.
#' @param ps.deriv.scale.rs For Weibull risk-set calibration:The derivative of \code{ps} with respect to the scale parameters of the Weibull risk-set calibration models.
# @param w A matrix of time points when measurements on the binary covariate were obtained.
# @param w.res A matrix of measurement results of the binary covariate. Each measurement corresponds to the time points in \code{w}
# @return Variance estimate for the log hazard-ratio of the binary exposure under  Weibull risk-set calibration model.
# @details DETAILS
# @examples 
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname CalcVar
#' @export 
CalcVarThetaWeibRS <- function(beta, etas.matrix, tm, event, ps.rs, ps.deriv.shape.rs, ps.deriv.scale.rs, w, w.res)
{
  n <- length(tm)
  b.vec <- Calcb(beta = beta, tm = tm, event = event, ps = ps.rs)
  nabla.etas.shape.Ubeta <- CalcUbetabeetaRS(beta = beta, tm = tm, event = event, ps = ps.rs, psDeriv = ps.deriv.shape.rs)
  nabla.etas.scale.Ubeta <- CalcUbetabeetaRS(beta = beta, tm = tm, event = event, ps = ps.rs, psDeriv = ps.deriv.scale.rs)
  nabla.etas.Ubeta <- c(rbind(nabla.etas.shape.Ubeta, nabla.etas.scale.Ubeta))/n

  hess.eta.inv <- ICweibHessSolvedRS(etas.matrix = etas.matrix, w = w, w.res = w.res, obs.tm = tm, event = event)
  grad.eta.pers <- ICweibGradRS(etas.matrix = etas.matrix, w = w, w.res = w.res,  obs.tm = tm, event = event)
   r.vec <- b.vec - nabla.etas.Ubeta%*%hess.eta.inv%*%t(grad.eta.pers)
  meat <- mean(r.vec^2)  # since beta is one-dimensional here 
  bread <- myFmyHess(beta, tm, event, ps.rs)/n
  var.beta <- (meat/(bread^2))/n
  return(var.beta)
}
# @title FUNCTION_TITLE
# @description FUNCTION_DESCRIPTION
# @param etas PARAM_DESCRIPTION
# @param w A matrix of time points when measurements on the binary covariate were obtained.
# @param w.res A matrix of measurement results of the binary covariate. Each measurement corresponds to the time points in \code{w}
# @return OUTPUT_DESCRIPTION
# @details DETAILS
# @examples 
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
# @seealso 
#  \code{\link[numDeriv]{hessian}}
# @rdname CalcVarEta
# @export 
# @importFrom numDeriv hessian
CalcVarEta <- function(etas,  w, w.res)
{
  n <- nrow(w)
  hess.etas.l.v <- (numDeriv::hessian(func = ICweibLik, x = etas, w = w, w.res = w.res))
  grad.eta.pers <- ICweibGrad(etas = etas, w = w, w.res = w.res)
  grad.eta <- 0
  for (j in 1:nrow(grad.eta.pers))
  {
    grad.eta <- grad.eta + grad.eta.pers[j,]%*%t(grad.eta.pers[j,])
  }    
  var.eta <-  solve(hess.etas.l.v)%*%(grad.eta)%*%solve(hess.etas.l.v)
    return(var.eta)
}
# @title Variance estimation and confidence intervals for the effect of the exposure under nonparametric calibration.
#' @description For nonparametric calibration, bootstrap calculations of the variance and confidence interval for the for the log hazard-ratio of the binary exposure.
# @param tm Vector of observed main event time or censoring time
# @param event Vector of censoring indicators. \code{1} for event \code{0} for censored
# @param w A matrix of time points when measurements on the binary covariate were obtained.
# @param w.res A matrix of measurement results of the binary covariate. Each measurement corresponds to the time points in \code{w}
#' @param BS For nonparametric calibration: Number of bootstrap iterations, Default: 100
#' @param CI For nonparametric calibration: Should the function return confidence intervals?, Default: T
#' @return For nonparametric calibration: Variance estimate and possibly confidence interval for the log hazard-ratio of the binary exposure under 
#' a nonparametric calibration model.
# @details DETAILS
# @examples 
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
# @seealso 
#  \code{\link[stats]{optimize}},\code{\link[stats]{cor}},\code{\link[stats]{quantile}}
#' @rdname CalcVar
#' @export 
#' @importFrom stats optimize var quantile
CalcVarNpmle <- function(tm, event, w, w.res, BS = 100, CI = T)
{
  n <- length(tm)
  beta.np.calib.bs <- vector(length = BS)
  for (j in 1:BS)
  {
    indices <- sample.int(n = n, size = n, replace = T)
    tm.bs <- tm[indices]
    event.bs <- event[indices]
    case.times.bs <- tm.bs[event.bs==1]
    w.bs <- w[indices,]
    w.res.bs <- w.res[indices,]
    fit.npmle.bs <- FitCalibNpmle(w = w.bs, w.res = w.res.bs)
    px.np.bs <- t(sapply(case.times.bs, CalcNpmleCalibP, w = w.bs, w.res =  w.res.bs, fit.npmle = fit.npmle.bs))
    px.np.bs[is.na(px.np.bs)] <- 0 # To avoid very rare errors
    beta.np.calib.bs[j] <- stats::optimize(f = CoxLogLikX,  tm = tm.bs, event = event.bs, ps = px.np.bs, 
                             interval = c(-50,50), maximum = T)$maximum
  }
  v.hat.npmle <- stats::var(beta.np.calib.bs[abs(beta.np.calib.bs) < 5])
  if (CI ==T)
  {
    ci.l <- stats::quantile(x = beta.np.calib.bs[abs(beta.np.calib.bs) < 5], probs = 0.025)
    ci.h <- stats::quantile(x = beta.np.calib.bs[abs(beta.np.calib.bs) < 5], probs = 0.975)
    return(list(v = v.hat.npmle, ci = c(ci.l, ci.h)))
  } else{
    return(list(v = v.hat.npmle))
}}
# @title Variance estimation and confidence intervals for the effect of the exposure under nonparametric  risk-set calibration models.
# @description Bootstrap calculations of the variance and confidence interval for the for the log hazard-ratio of the binary exposure
# under nonparametric risk-set calibration models.
# @param tm Vector of observed main event time or censoring time
# @param event Vector of censoring indicators. \code{1} for event \code{0} for censored
# @param w A matrix of time points when measurements on the binary covariate were obtained.
# @param w.res A matrix of measurement results of the binary covariate. Each measurement corresponds to the time points in \code{w}
# @param BS Number of bootstrap iterations, Default: 100
# @param CI Should the function return confidence intervals?, Default: T
# @return Variance estimate and possibly confidence interval for the log hazard-ratio of the binary exposure under 
# nonparametric risk-set calibration models.
# @details DETAILS
# @examples 
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
# @seealso 
#  \code{\link[stats]{optimize}},\code{\link[stats]{cor}},\code{\link[stats]{quantile}}
#' @rdname CalcVar
#' @export 
#' @importFrom stats optimize var quantile
CalcVarNpmleRS <- function(tm, event, w, w.res, BS = 100, CI =T)
{
  n <- length(tm)
  beta.np.calib.rs.bs <- vector(length = BS)
  for (j in 1:BS)
  {
    #  cat("j = ", j)
    indices <- sample.int(n = n, size = n, replace = T)
    tm.bs <- tm[indices]
    event.bs <- event[indices]
    case.times.bs <- tm.bs[event.bs==1]
    w.bs <- w[indices,]
    w.res.bs <- w.res[indices,]
    px.np.rs.bs <- t(sapply(case.times.bs, CalcNpmleRSP, w = w.bs, w.res =  w.res.bs, obs.tm = tm.bs))
    px.np.rs.bs[is.na(px.np.rs.bs)] <- 0 # To avoid very rare errors
    beta.np.calib.rs.bs[j] <- stats::optimize(f = CoxLogLikX,  tm = tm.bs, event = event.bs, ps = px.np.rs.bs, 
                                    interval = c(-50,50), maximum = T)$maximum
  }
  v.hat.npmle.rs <- stats::var(beta.np.calib.rs.bs[abs(beta.np.calib.rs.bs) < 5])
  if (CI ==T)
  {
    ci.l <- stats::quantile(x = beta.np.calib.rs.bs[abs(beta.np.calib.rs.bs) < 5], probs = 0.025)
    ci.h <- stats::quantile(x = beta.np.calib.rs.bs[abs(beta.np.calib.rs.bs) < 5], probs = 0.975)
    return(list(v = v.hat.npmle.rs, ci = c(ci.l, ci.h)))
  } else{
    return(list(v = v.hat.npmle.rs))
}}

# CalcVarThetaCox <- function(beta, etas, tm, event, ps, ps.deriv.shape, ps.deriv.scale, w, w.res)
# {
#   n <- length(tm)
#   b.vec <- Calcb(beta = beta, tm = tm, event = event, ps = ps)
#   nabla.eta.shape.Ubeta <- CalcUbetabeeta(beta = beta, tm = tm, event = event, ps = ps, psDeriv = ps.deriv.shape)
#   nabla.eta.scale.Ubeta <- CalcUbetabeeta(beta = beta, tm = tm, event = event, ps = ps, psDeriv = ps.deriv.scale)
#   nabla.eta.Ubeta <- c(nabla.eta.shape.Ubeta, nabla.eta.scale.Ubeta)/n
#   hess.etas.l.v <- (hessian(func = ICweibLik, x = etas, w = w, w.res = w.res))
#   grad.eta.pers <- ICweibGrad(etas = etas, w = w, w.res = w.res)
#   r.vec <- b.vec - nabla.eta.Ubeta%*%solve(hess.etas.l.v)%*%t(grad.eta.pers)
#   meat <- mean(r.vec^2)  # since beta is one-dimensional here 
#   bread <- myFmyHess(beta, tm, event, ps)/n
#   var.beta <- (meat/(bread^2))/n
#   return(var.beta)
# }
# 
# CalcVarThetaCoxRS <- function(beta, etas.matrix, tm, event, ps.rs, ps.deriv.shape.rs, ps.deriv.scale.rs, w, w.res)
# {
#   n <- length(tm)
#   b.vec <- Calcb(beta = beta, tm = tm, event = event, ps = ps.rs)
#   nabla.etas.shape.Ubeta <- CalcUbetabeetaRS(beta = beta, tm = tm, event = event, ps = ps.rs, psDeriv = ps.deriv.shape.rs)
#   nabla.etas.scale.Ubeta <- CalcUbetabeetaRS(beta = beta, tm = tm, event = event, ps = ps.rs, psDeriv = ps.deriv.scale.rs)
#   nabla.etas.Ubeta <- c(rbind(nabla.etas.shape.Ubeta, nabla.etas.scale.Ubeta))/n
#   
#   hess.eta.inv <- ICweibHessSolvedRS(etas.matrix = etas.matrix, w = w, w.res = w.res, tm = tm, event = event)
#   grad.eta.pers <- ICweibGradRS(etas = etas.matrix, w = w, w.res = w.res,  tm = tm, event = event)
#   r.vec <- b.vec - nabla.etas.Ubeta%*%hess.eta.inv%*%t(grad.eta.pers)
#   meat <- mean(r.vec^2)  # since beta is one-dimensional here 
#   bread <- myFmyHess(beta, tm, event, ps.rs)/n
#   var.beta <- (meat/(bread^2))/n
#   return(var.beta)
# }
