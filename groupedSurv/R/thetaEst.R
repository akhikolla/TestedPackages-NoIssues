# function to estimate thetas for the covarate

.grad_NF <- function(Params, Xmatrix, Kivec, Deltavec, ntps) {
  .Call("_groupedSurv_grad_NF", PACKAGE = "groupedSurv", Params, Xmatrix, Kivec, 
    Deltavec, ntps)
}

.logLike_NF <- function(Params, Xmatrix, Kivec, Deltavec, ntps) {
  .Call("_groupedSurv_logLike_NF", PACKAGE = "groupedSurv", Params, Xmatrix, Kivec, 
    Deltavec, ntps)
}

thetaEst <- function(Z=NULL, gtime, delta, method="BFGS")
{
  if(sum(is.infinite(gtime)) >= 1)
     ntps <- nlevels(as.factor(gtime)) - 1
  else
     ntps <- nlevels(as.factor(gtime))

	alphaIG <- runif(ntps, 0, 1)
  # alphaIG <- alphaEstFam(Dtime, Event) +0.0000000000000000005
  thetaIG <- runif(ncol(Z), 0, 1)
  # cat('thetaTG: ', thetaIG, '\n')
  ThetaIG <- c(alphaIG, thetaIG)

  Z <- as.matrix(Z,ncol=ncol(Z))
  Est <- optim(par = ThetaIG, fn = .logLike_NF, gr = .grad_NF, Xmatrix = Z, Kivec = gtime, 
    Deltavec = delta, ntps = ntps, method = method, control = list(fnscale = -1))$par
  ## testing other methods for optim Est <- optim(par = ThetaIG, fn = .logLike_NF,
  ## gr = .grad_NF, Xmatrix = Z, Kivec = Dtime, Deltavec = Event, ntps = ntps,
  ## method = 'CG', control = list(fnscale = -1))$par
  
  thetaest <- NULL
  thetaest$alpha <- exp(-exp(Est[1:ntps]))
  thetaest$theta <- Est[(ntps + 1):length(Est)]
  
  thetaest
}

















