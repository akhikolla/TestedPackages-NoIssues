# Weibull auxiliary functions for likelihood, gradient and hessian calculations for Weibull fitting from interval-censored data
ICweibLik <- function(etas,w,w.res)
{
  lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
  L <- lr.for.fit[,1]
  R <- lr.for.fit[,2]
  eta.shape <- etas[1]
  eta.scale <- etas[2]
  pers <- log(exp(-(L/eta.scale)^eta.shape)-exp(-(R/eta.scale)^eta.shape))
  s <- sum(pers)
  return(s)
}

ICweibGrad <- function(etas,w,w.res)
{
  lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
  L <- lr.for.fit[,1]
  R <- lr.for.fit[,2]
  LRwithin <- (R <Inf & L>0)
  Lzero <- (L==0 & R<Inf)
  Rinf <- (R==Inf & L >0)
  LRout <- (R==Inf & L ==0)
  eta.shape <- etas[1]
  eta.scale <- etas[2]
  pers.lik <- (exp(-(L/eta.scale)^eta.shape)-exp(-(R/eta.scale)^eta.shape))
  deriv.eta.scale.pers <- deriv.eta.shape.pers <- vector(length = length(L))
  deriv.eta.scale.pers[LRwithin] <- (eta.shape/(pers.lik[LRwithin]*(eta.scale^(eta.shape+1)))) * ((L[LRwithin]^eta.shape)*exp(-(L[LRwithin]/eta.scale)^eta.shape) - (R[LRwithin]^eta.shape)*exp(-(R[LRwithin]/eta.scale)^eta.shape))
  deriv.eta.scale.pers[Rinf] <- (eta.shape/(pers.lik[Rinf]*eta.scale^(eta.shape+1))) * ((L[Rinf]^eta.shape)*exp(-(L[Rinf]/eta.scale)^eta.shape))
  deriv.eta.scale.pers[Lzero] <- (eta.shape/(pers.lik[Lzero]*eta.scale^(eta.shape+1))) * (- (R[Lzero]^eta.shape)*exp(-(R[Lzero]/eta.scale)^eta.shape))
  deriv.eta.scale.pers[LRout] <- 0
  deriv.eta.shape.pers[LRwithin] <- (1/(pers.lik[LRwithin]*eta.scale^eta.shape)) * (log(R[LRwithin]/eta.scale)*exp(-(R[LRwithin]/eta.scale)^eta.shape)*(R[LRwithin]^eta.shape) - log(L[LRwithin]/eta.scale)*exp(-(L[LRwithin]/eta.scale)^eta.shape)*(L[LRwithin]^eta.shape))
  deriv.eta.shape.pers[Rinf] <- (1/(pers.lik[Rinf] *eta.scale^eta.shape)) * (- log(L[Rinf] /eta.scale)*exp(-(L[Rinf]/eta.scale)^eta.shape)*(L[Rinf] ^eta.shape))
  deriv.eta.shape.pers[Lzero] <- (1/(pers.lik[Lzero]*eta.scale^eta.shape)) * (log(R[Lzero]/eta.scale)*exp(-(R[Lzero]/eta.scale)^eta.shape)*(R[Lzero]^eta.shape))
  deriv.eta.shape.pers[LRout] <- 0
  return(cbind(deriv.eta.shape.pers,deriv.eta.scale.pers))
}


ICweibLikRS <- function(etas.matrix.in.vector, w, w.res, obs.tm, event)
{
  r <- sum(event)
  s <- 0
  if(r != length(etas.matrix.in.vector)/2) stop('Someting is wrong!')
  etas.matrix <- matrix(nrow = r, ncol = 2, etas.matrix.in.vector, byrow = F)
  event.index <- which(event==1)
  lr.for.fit.all <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
  u.risk <- vector(length = r)
  for(j in 1:r)
  {
    point <- obs.tm[event.index[j]]
        lr.for.fit <- lr.for.fit.all[obs.tm>=point,]  # Keep only observations in the risk set
  L <- lr.for.fit[,1]
  R <- lr.for.fit[,2]
  eta.shape <- etas.matrix[j,1]
  eta.scale <- etas.matrix[j,2]
  u.risk[j] <- sum(log(exp(-(L/eta.scale)^eta.shape)-exp(-(R/eta.scale)^eta.shape)))
  }
  s <- sum(u.risk)
  return(s)
}


ICweibGradRS <- function(etas.matrix, w, w.res, obs.tm, event)
{
  n <- length(event)
  r <- sum(event)
  if(r != nrow(etas.matrix)) stop('Someting is wrong!')
  mat.shape.back <- matrix(nrow = n, ncol = r)
  mat.scale.back <- matrix(nrow = n, ncol = r)
  event.index <- which(event==1)
  lr.for.fit <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
  L <- lr.for.fit[,1]
  R <- lr.for.fit[,2]
  LRwithin <- (R <Inf & L>0)
  Lzero <- (L==0 & R<Inf)
  Rinf <- (R==Inf & L >0)
  LRout <- (R==Inf & L ==0)
  for(j in 1:r)
  {
  point <- obs.tm[event.index[j]]
  eta.shape <- etas.matrix[j, 1]
  eta.scale <- etas.matrix[j, 2]
  in.risk.set <- obs.tm>=point
  pers.lik <- (exp(-(L/eta.scale)^eta.shape)-exp(-(R/eta.scale)^eta.shape))
  mat.scale.back[LRwithin, j] <- (eta.shape/(pers.lik[LRwithin]*(eta.scale^(eta.shape+1)))) * ((L[LRwithin]^eta.shape)*exp(-(L[LRwithin]/eta.scale)^eta.shape) - (R[LRwithin]^eta.shape)*exp(-(R[LRwithin]/eta.scale)^eta.shape))
  mat.scale.back[Rinf, j] <- (eta.shape/(pers.lik[Rinf]*eta.scale^(eta.shape+1))) * ((L[Rinf]^eta.shape)*exp(-(L[Rinf]/eta.scale)^eta.shape))
  mat.scale.back[Lzero, j] <- (eta.shape/(pers.lik[Lzero]*eta.scale^(eta.shape+1))) * (- (R[Lzero]^eta.shape)*exp(-(R[Lzero]/eta.scale)^eta.shape))
  mat.scale.back[LRout | !in.risk.set, j] <- 0
  mat.shape.back[LRwithin, j] <- (1/(pers.lik[LRwithin]*eta.scale^eta.shape)) * (log(R[LRwithin]/eta.scale)*exp(-(R[LRwithin]/eta.scale)^eta.shape)*(R[LRwithin]^eta.shape) - log(L[LRwithin]/eta.scale)*exp(-(L[LRwithin]/eta.scale)^eta.shape)*(L[LRwithin]^eta.shape))
  mat.shape.back[Rinf, j] <- (1/(pers.lik[Rinf] *eta.scale^eta.shape)) * (- log(L[Rinf] /eta.scale)*exp(-(L[Rinf]/eta.scale)^eta.shape)*(L[Rinf] ^eta.shape))
  mat.shape.back[Lzero, j] <- (1/(pers.lik[Lzero]*eta.scale^eta.shape)) * (log(R[Lzero]/eta.scale)*exp(-(R[Lzero]/eta.scale)^eta.shape)*(R[Lzero]^eta.shape))
  mat.shape.back[LRout | !in.risk.set, j] <- 0
  }
  mat.grad.back <- matrix(nrow = n, ncol = 2*r)
  mat.grad.back[,seq(1, 2*r, 2)] <- mat.shape.back
  mat.grad.back[,seq(2, 2*r, 2)] <- mat.scale.back
  return(mat.grad.back)
}

ICweibLikLR <- function(weib.params, L,R)
{
  eta.shape <- weib.params[1]
  eta.scale <- weib.params[2]
  u.risk <- sum(log(exp(-(L/eta.scale)^eta.shape)-exp(-(R/eta.scale)^eta.shape)))
  return(u.risk)
}

# @importFrom numDeriv hessian
ICweibHessSolvedRS <- function(etas.matrix, w, w.res, obs.tm, event)
{
  n <- length(event)
  r <- sum(event)
  event.index <- which(event==1)
  lr.for.fit.all <- as.data.frame(FindIntervalCalibCPP(w = w, wres = w.res))
  hess.etas.solved <- matrix(nrow = 2*r, ncol = 2*r, 0)
for (j in 1:r)
{
  weib.param <- etas.matrix[j,]
  point <- obs.tm[event.index[j]]
  lr.for.fit <- lr.for.fit.all[obs.tm>=point,]  # Keep only observations in the risk set
  L <- lr.for.fit[,1]
  R <- lr.for.fit[,2]
  hess.etas.solved[(2*j-1):(2*j),(2*j-1):(2*j)] <- solve(numDeriv::hessian(func = ICweibLikLR, x = weib.param, L = L, R =R))
}
 return(hess.etas.solved) 
}
