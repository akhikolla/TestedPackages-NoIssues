



GM <- function(
  X,
  Xw = rep(1.0, ncol(X)),
  alpha = numeric(0),
  mu = matrix(ncol = 0, nrow = 0),
  sigma = matrix(ncol = 0, nrow = 0),
  G = 5L,
  convergenceEPS = 1e-5,
  alphaEPS = 0,
  eigenRatioLim = Inf,
  maxIter = 1000L,
  maxCore = 7L,
  tlimit = 3600,
  verbose = TRUE,
  updateAlpha = TRUE,
  updateMean = TRUE,
  updateSigma = TRUE,
  paraConvergeMaxErr = FALSE,
  loglikehoodConverge = FALSE,
  loglikehoodConvergeBlock = 10
  )
{
  if(!is.finite(eigenRatioLim)) eigenRatioLim = 0;
  rst = paraGmm(
    X,
    Xw,
    G,
    alpha,
    mu,
    sigma,
    eigenRatioLim,
    convergenceEPS,
    alphaEPS,
    maxIter,
    tlimit,
    verbose,
    maxCore,
    updateAlpha,
    updateMean,
    updateSigma,
    paraConvergeMaxErr,
    loglikehoodConverge,
    loglikehoodConvergeBlock
  )
  rst$clusterMember = aggregate(list(1L : ncol(X)), list(rst$clusterMember), function(x) x)[[2]]
  rst
}




GMcw <- function(
  X,
  Xw = rep(1.0, ncol(X)),
  alpha = numeric(0),
  mu = matrix(ncol = 0, nrow = 0),
  sigma = matrix(ncol = 0, nrow = 0),
  G = 5L,
  convergenceEPS = 1e-5,
  alphaEPS = 0,
  eigenRatioLim = Inf,
  maxIter = 1000L,
  maxCore = 7L,
  tlimit = 3600,
  verbose = TRUE)
{
  if(!is.finite(eigenRatioLim)) eigenRatioLim = 0;
  rst = paraGmmCW(
    X,
    Xw,
    G,
    alpha,
    mu,
    sigma,
    eigenRatioLim,
    convergenceEPS,
    alphaEPS,
    maxIter,
    tlimit,
    verbose,
    maxCore
  )
  rst$clusterMember = aggregate(list(1L : ncol(X)), list(rst$clusterMember), function(x) x)[[2]]
  rst
}




GMfj <- function(
  X,
  Xw = rep(1.0, ncol(X)),
  alpha = numeric(0),
  mu = matrix(ncol = 0, nrow = 0),
  sigma = matrix(ncol = 0, nrow = 0),
  G = 5L,
  Gmin = 2L,
  convergenceEPS = 1e-5,
  alphaEPS = 0,
  eigenRatioLim = Inf,
  maxIter = 1000L,
  maxCore = 7L,
  tlimit = 3600,
  verbose = TRUE)
{
  if(!is.finite(eigenRatioLim)) eigenRatioLim = 0;
  rst = paraGmmFJ(
    X,
    Xw,
    G,
    Gmin,
    alpha,
    mu,
    sigma,
    eigenRatioLim,
    convergenceEPS,
    alphaEPS,
    maxIter,
    tlimit,
    verbose,
    maxCore)
  # print(rst$clusterMember)
  rst$clusterMember = aggregate(list(1L : ncol(X)), list(rst$clusterMember), function(x) x)[[2]]
  rst
}








































































