#### CalcWeibullCalibP function
### NEED TO WRITE EXPLANATION
## The function takes the following
## w - a matrix. Each row is observation and each column is questionnaire time in the interval. w equal to Inf once
# an observation is censore/had the event
## w.res - a matrix of the same dimensions as w. Equal to the x(t) at time w. For example second column is 
# second questionnaire result for all participents.
## point - scalar. The time of the risk set in the main analysis. In terms of the paper, t.
## weib.params - the shape and scale parameters from the Weibull calibration fitting to the interval-cenosed time to exposure/treatment
###
# The function returns a vector with individual predictions for P(X(t)=1|history(time t)). 
# For observations with X(a(t))=1 the above probability is 1 by definition and this is what the
# function returns for them.
# The function calculates prediction for all observations, even though predictions for observations outside 
# the risk set are not used
#### The following functions is used: CalcAuxatPoint (R function)
# @importFrom stats pweibull
CalcWeibullCalibPderivScale <- function(w, w.res, point, weib.params)
{
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  weib.shape <- weib.params[1]
  weib.scale <- weib.params[2]
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  surv.at.point <- stats::pweibull(point, shape = weib.shape,scale = weib.scale, lower.tail = F)
  surv.at.a.point <- stats::pweibull(a.point[p.point==0], shape = weib.shape, scale = weib.scale, lower.tail = F)
  deriv.scale <- vector(length=length(p.point))
  deriv.scale[p.point==0] <- -(weib.shape/(weib.scale^(weib.shape + 1))) * (point^weib.shape - a.point[p.point==0]^weib.shape) * (surv.at.point/surv.at.a.point)
  deriv.scale[p.point>0] <- 0
  return(deriv.scale)
}

# @importFrom stats pweibull
CalcWeibullCalibPderivShape <- function(w, w.res, point, weib.params)
{
  weib.shape <- weib.params[1]
  weib.scale <- weib.params[2]
  lr.for.lik <- CalcAuxAtPoint(w,w.res,point = point)
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  surv.at.point <- stats::pweibull(point, shape = weib.shape,scale = weib.scale, lower.tail = F)
  surv.at.a.point <- stats::pweibull(a.point[p.point==0], shape = weib.shape, scale = weib.scale, lower.tail = F)
  deriv.shape <- vector(length=length(p.point))
  deriv.shape[p.point==0] <- (1/(weib.scale^weib.shape)) * (log(a.point[p.point==0]/weib.scale)*a.point[p.point==0]^weib.shape - 
                                                             log(point/weib.scale)*point^weib.shape ) * (surv.at.point/surv.at.a.point)
  deriv.shape[p.point>0] <- 0
  deriv.shape[a.point==0] <- -(1/(weib.scale^weib.shape)) * log(point/weib.scale)*point^weib.shape * surv.at.point
  return(deriv.shape)
}


#########################################################################################################################
# @importFrom stats pweibull
CalcWeibullCalibPderivShapeRS <- function(w, w.res, obs.tm, event, weib.rs.params)
{
  r <- sum(event)
  n <- length(event)
  event.index <- which(event==1)
  deriv.shape <- matrix(nrow = n, ncol = r)
  for (j in 1:r)
  {
  point <- obs.tm[event.index[j]]
  weib.rs.shape <- weib.rs.params[j, 1]
  weib.rs.scale <- weib.rs.params[j, 2]
  in.risk.set <- obs.tm>=point
  lr.for.lik <- CalcAuxAtPoint(w, w.res, point = point)
  a.point <- lr.for.lik$a.point
  p.point <- lr.for.lik$x.one
  surv.at.point <- stats::pweibull(point, shape = weib.rs.shape, scale = weib.rs.scale, lower.tail = F)
  surv.at.a.point <- stats::pweibull(a.point[p.point==0 & in.risk.set], shape = weib.rs.shape, scale = weib.rs.scale, lower.tail = F)
  deriv.shape[p.point==0 & in.risk.set, j] <- (1/(weib.rs.scale^weib.rs.shape)) * (log(a.point[p.point==0 & in.risk.set]/weib.rs.scale) * 
                                                a.point[p.point==0 & in.risk.set]^weib.rs.shape - 
                                                  log(point/weib.rs.scale)*point^weib.rs.shape ) *
                                                               (surv.at.point/surv.at.a.point)
  deriv.shape[p.point>0, j] <- 0
  deriv.shape[a.point==0 & in.risk.set, j] <- -(1/(weib.rs.scale^weib.rs.shape)) * log(point/weib.rs.scale)*point^weib.rs.shape * surv.at.point
  deriv.shape[!in.risk.set, j] <- 0
  }
  return(t(deriv.shape))
}

#########################################################################################################################

# @importFrom stats pweibull
CalcWeibullCalibPderivScaleRS <- function(w, w.res, obs.tm, event, weib.rs.params)
{
  r <- sum(event)
  n <- length(event)
  event.index <- which(event==1)
  deriv.scale <- matrix(nrow = n, ncol = r)
  for (j in 1:r)
  {
    point <- obs.tm[event.index[j]]
    weib.rs.shape <- weib.rs.params[j, 1]
    weib.rs.scale <- weib.rs.params[j, 2]
    in.risk.set <- obs.tm>=point
    lr.for.lik <- CalcAuxAtPoint(w, w.res, point = point)
    a.point <- lr.for.lik$a.point
    p.point <- lr.for.lik$x.one
    surv.at.point <- stats::pweibull(point, shape = weib.rs.shape,scale = weib.rs.scale, lower.tail = F)
    surv.at.a.point <- stats::pweibull(a.point[p.point==0 & in.risk.set], shape = weib.rs.shape, scale = weib.rs.scale, lower.tail = F)
    deriv.scale[p.point==0 & in.risk.set, j] <- -(weib.rs.shape/(weib.rs.scale^(weib.rs.shape + 1))) * 
    (point^weib.rs.shape - a.point[p.point==0 & in.risk.set]^weib.rs.shape) * (surv.at.point/surv.at.a.point)
    deriv.scale[p.point>0 , j] <- 0
    deriv.scale[!in.risk.set, j] <- 0
  }
  return(t(deriv.scale))
}
