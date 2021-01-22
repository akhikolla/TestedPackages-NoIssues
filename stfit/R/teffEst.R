## Copyright 2011 – 2020  Center for Survey Statistics and Methodology at Iowa State University   All Rights Reserved
#' STFIT Temporal Effect Estimation
#'
#' @param ids ids for 'group', for data with repeated measurement over years, year is ids; 
#' for pixels belong to certain clusters, cluster is ids.
#' @param doy vecotr of DOY (day of the year)
#' @param rmat residual matrix with rows corresponding to \code{doy} and columns corresponding to pixel index
#' @param h.cov bandwidth for temporal covariance estimation; ignored if \code{weight.cov} is supplied
#' @param h.sigma2 bandwidth for temporal variance estimation
#' @param weight.cov weight vector for temporal covariance estimation
#' @param weight.sigma2 weight vector for temporal variance estimation
#' @param pve percentage of variance explained; used for number of eigen values selection. Default is 0.99.
#' @param doyeval a vector of DOY on which to get the temporal imputation
#' @param t.grid a vector of grid points on which to calculate the temporal covariance function
#' @param t.grid.num number of grid points to use for temporal covariance estimation. Ignored if \code{t.grid} is given.
#' @param var.est Whether to estimate the variance of the temporal effect. Default is FALSE.
#'
#' @return List of length 2 with entries:
#'   \itemize{
#'     \item teff_array: 3-d array with first dimention 'ids', second dimention 'doy' and third
#'     dimention pixel index.
#'     \item teff_var_array: same structure as \code{teff_array} if \code{var.est} is TRUE,
#'     otherwise NULL.
#'   }
#' @export
#' @example 
#' \donttest{
#' library(doParallel)
#' library(raster)
#' library(rasterVis)
#' library(RColorBrewer)
#' dfB = landsat106[landsat106$year >= 2000,]
#' matB = as.matrix(dfB[,-c(1:2)])
#' year = dfB$year
#' doy = dfB$doy
#' ## Speed up the calculation by using multi-cores if doParallel package is installed
#' if(require(doParallel))
#'   registerDoParallel(1)
#' meanest = meanEst(doy, matB, 1:365)
#' rmat = matB - meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]
#' teffres = teffEst(year, doy, rmat, doyeval = meanest$doyeval)
#' plot(teffres$teff_array[1,,1],type = 'l')
#' }
teffEst <- function(ids, doy, rmat,
                    doyeval = seq(min(doy), max(doy)), h.cov = 100, h.sigma2 = 300,
                    weight.cov = NULL, weight.sigma2 = NULL,
                    pve = 0.99, t.grid = NULL, t.grid.num = 50, var.est = FALSE){
  message("Estimating the temporal effect...\n")
  if(is.null(weight.cov))
    weight.cov = weightVector(h.cov)
  if(is.null(weight.sigma2))
    weight.sigma2 = weightVector(h.sigma2)
  yeareval = sort(unique(ids))
  if(is.null(t.grid)){
    if(t.grid.num <= length(doyeval))
      t.grid = doyeval[unique(round(seq(1, length(doyeval), length.out=t.grid.num)))] else
        stop("t.grid.num is bigger than length of doyeval.")
  }
  acomb <- function(...) abind::abind(..., along=3)
  i = NULL
  teffres = foreach(i = 1:ncol(rmat)) %dopar% {
    resid = rmat[,i]
    nnaidx = !is.na(resid)
    if(sum(nnaidx) == 0)
      return(matrix(0, length(yeareval), length(doyeval)))
    R0.hat = lc_cov_1d_est(ids[nnaidx], doy[nnaidx], resid[nnaidx], weight.cov, t.grid)
    sigma2 = llreg(doy[nnaidx], resid[nnaidx]^2, x.eval = t.grid, h = h.sigma2)
    
    ## eigen decomposition
    ev = eigen(R0.hat)
    ev$values[ev$values < 0] = 0
    ev.idx = max(which.min(cumsum(ev$values)/sum(ev$values) < pve), 2) ## select at least 2 eigen components
    ## message("The first ", ev.idx, " eigen values are used...\n")
    ev.vec = ev$vectors[, 1:ev.idx, drop=FALSE]
    ev.val = ev$values[1:ev.idx]
    ev.vec = phiInterp(doyeval, ev.vec, t.grid)
    
    idx1 = round(length(t.grid)/4):round(length(t.grid)/4*3)
    nugg = max(mean((sigma2[idx1]-diag(R0.hat)[idx1])), 1e-06)
    if(!all(yeareval %in% ids[nnaidx])){
      res = PACE1d(ids[nnaidx], doy[nnaidx], resid[nnaidx], ev.vec, nugg, ev.val, doyeval, var.est = var.est)
      tmat = matrix(0, length(yeareval), length(doyeval))
      tmat[which(yeareval %in% ids[nnaidx]),] = res$tmat
      if(var.est){
        tvarmat = matrix(0, length(yeareval), length(doyeval))
        tvarmat[which(yeareval %in% ids[nnaidx]),] = res$tvarmat
      } else {
        tvarmat = NULL
      }
    } else{
      res = PACE1d(ids[nnaidx], doy[nnaidx], resid[nnaidx], ev.vec, nugg, ev.val, doyeval, yeareval, var.est = var.est)
      tmat = res$tmat
      if(var.est){
        tvarmat = res$tvarmat
      } else {
        tvarmat = NULL
      }
    }
    return(list(tmat = tmat, tvarmat = tvarmat))
  }
  teff_array = abind::abind(lapply(teffres, function(x) x[[1]]), along = 3)
  dimnames(teff_array) = list(yeareval, doyeval, 1:ncol(rmat))
  teff_var_array = NULL
  if(var.est){
    teff_var_array = abind::abind(lapply(teffres, function(x) x[[2]]), along = 3)
    dimnames(teff_var_array) = list(yeareval, doyeval, 1:ncol(rmat))
  }
  return(list(teff_array = teff_array, teff_var_array = teff_var_array))
}
