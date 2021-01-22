## Copyright 2011 – 2020  Center for Survey Statistics and Methodology at Iowa State University   All Rights Reserved
#' STFIT Mean Estimation
#' 
#' The function is used for pixel-wise mean estimation.
#' 
#' There are several predefined methods for mean estimation: \code{smooth_spline}, 
#' \code{llreg}, \code{lpreg} and \code{spreg}. User can use \code{opt$get()} to check
#' the current registered method and use \code{opt$set()} function to set the method.
#' For exmaple, one can run \code{opt$set(smooth_spline)} first and then run the 
#' \code{meanEst} function to use smoothing spline regression for mean eatimation.
#' User can also customize the methods for mean estimation. For example, mean estimation
#' through fourier basis expansion:
#' 
#' \preformatted{
#' .X = fda::eval.basis(1:365, fda::create.fourier.basis(rangeval=c(0,365), nbasis=11))
#' customfun <- function(x, y, x.eval=1:365, minimum.num.obs = 10){
#'   nonna.idx = !is.na(y)
#'   if(sum(nonna.idx) < minimum.num.obs)
#'     return(rep(NA, 365))
#'   ## lmfit = lm.fit(.X[unlist(lapply(x, function(x) which(x == x.eval))),], y[nonna.idx])
#'   lmfit = lm.fit(.X[x[nonna.idx],], y[nonna.idx])
#'   return(.X[x.eval,] %*% lmfit$coefficient)
#' }
#' stfit::opts_stfit$set(temporal_mean_est = customfun)
#' }
#' 
#' @param doy vector of day of year (DOY) index
#' @param mat data matrix. Each row contains a row stacked image pixel values.
#' @param doyeval a vector of DOY on which to get the mean imputation
#' @param outlier.tol the tolerance value in defining an image as outlier. The percent of 
#' outlier pixels in an image exceed this value is regarded as outlier image which will not
#' be used in temporal mean estimation.
#' @param msk an optional logistic vector. TRUE represent the corresponding pixel is always missing.
#' @param cluster an optional vector defining clusters of pixels. If NULL, mean estimation
#' is conducted on each pixel, otherwise all pixels from the same cluster are combined for
#' mean estimation.
#' @param minimum.num.obs minimum number of observations needed for mean estimation. Too few observations
#' may lead to big estimation error.
#' @param redo whether to recalculate the mean estimation if there is an outlier (only redo once).
#' @param clipRange vector of length 2, specifying the minimum and maximum values of the prediction value
#' @param clipMethod "nnr" or "truncate". "nnr" uses average of nearest neighbor pixels to impute;
#' "truncate use the clipRange value to truncate.
#' @param img.nrow number of rows for an image, only used when 'clipMethod' is "nnr"
#' @param img.ncol number of columns for an image, only used when 'clipMethod' is "nnr"
#' 
#' @return a list containing the following entries:
#'   \itemize{
#'     \item doyeval: same as input \code{doyeval}
#'     \item meanmat: estimated mean matrix, with number of rows equals length of \code{doyeval}
#'      and number of columns equal \code{ncol(mat)}
#'      \item idx: a list of image indexes 
#'        \itemize{
#'        \item idx.allmissing: completely missing image indexes,
#'        \item idx.partialmissing: partially observed image indexes,
#'        \item idx.fullyobserved: fully observed image indexes,
#'        \item idx.outlier: outlier image indexes.
#'      }
#'     \item outlier: a list of image outliers information
#'       \itemize{
#'         \item outidx: index of the outlier image
#'         \item outpct: percentage of outlier pixels corresponding to \code{outidx},
#'         \item outlst: a list of the same length as \code{outidx}, with each list the missing pixel index.
#'      }
#'   }
#' 
#' @export
#' 
#' @example
#' \donttest{
#' dfB = landsat106[landsat106$year >= 2000,]
#' matB = as.matrix(dfB[,-c(1:2)])
#' year = dfB$year
#' doy = dfB$doy
#' ## Speed up the calculation by using multi-cores if doParallel package is installed
#' if(require(doParallel))
#'   registerDoParallel(8)
#' meanest = meanEst(doy, matB, 1:365)
#' idx = round(c(seq(1, 365, length.out = 16)))
#' landsatVis(meanest$meanmat[idx,], names.attr = paste0("DOY", idx))
#' }
meanEst <- function(doy, mat,
                    doyeval = seq(min(doy), max(doy)), 
                    msk = rep(FALSE, ncol(mat)), outlier.tol = 0.5, minimum.num.obs = 4,
                    cluster = NULL, redo = TRUE, clipRange = c(-Inf, Inf), clipMethod = c("truncate", "nnr"),
                    img.nrow=NULL, img.ncol=NULL){
  clipMethod = match.arg(clipMethod)
  idx = 1:length(doy) ## idx is the index of image, 1, 2, 3,...
  N = ncol(mat) ## number of pixels
  M = length(doyeval) ## number of doy to evaluate on
  if(!all(doy %in% doyeval)){
    stop("'doyeval' should be a superset of 'doy")
  }
  if(length(doy) != nrow(mat))
    stop("doy and mat dimension do not match!")
  if(length(msk) != N)
    stop("msk and mat dimension do not match!")
  
  ################################################################
  ######### Cluster/Pixel-wise overall mean estimation ###########
  ################################################################
  ## mean.mat: columns are pixel index, rows are doy index (ex. 365 x 961)
  mean.mat = matrix(NA, M, N)
  if(is.null(cluster)){
    mean.mat[, !msk] = .pixelwiseMeanEst(doy, mat[, !msk], doyeval, minimum.num.obs, NULL)
  } else {
    mean.mat[, !msk] = .pixelwiseMeanEst(doy, mat[, !msk], doyeval, minimum.num.obs, cluster[!msk])
  }
  
  
  ## there will be a problem if some pixel has very few observation
  ## so check whether the missing pattern matches with 'msk'
  if(!all(msk == getMask(mean.mat)))
    warning(paste0("Some pixels which are not covered by mask have less than ", 
                   minimum.num.obs, " observations!"))
  
  ###########################################################
  ######### Outlier detection using the residuals ###########
  ###########################################################
  ## residual matrix
  resid.mat = mat - mean.mat[unlist(lapply(doy, function(x,y) which(y == x), y = doyeval)),]
  outlier.res = outlier(resid.mat)
  oidx = outlier.res$outpct > outlier.tol
  ## outlier image indexes
  idx0 = outlier.res$outidx[oidx]
  
  pct_missing = apply(mat, 1, function(x) {sum(is.na(x))/N})
  ## Index
  idx1 = idx[pct_missing == 1] ## all missing images indexes;
  ## partially missing image indexes after removing outliers;
  idx2 = setdiff(idx[(pct_missing > 0) & (pct_missing < 1)], idx0)
  ## no missing images indexes after removing outliers;
  idx3 = setdiff(idx[pct_missing == 0], idx0)
  
  ## Outliers are relaced with NA
  message("Outlier pixels are replaced with NAs...\n")
  ## use NA instead of the outlier pixel values both in original data
  for(i in 1:length(outlier.res$outidx)){
    mat[outlier.res$outidx[i], outlier.res$outlst[[i]]] = NA
  }
  
  ##########################################
  ######### Redo mean estimation ###########
  ##########################################
  ## Redo pixel-wise temporal trend estimation after removing outliers
  if (length(idx0) > 0)
    message(paste0(
      length(idx0),
      " outlier images are removed in the mean estimation procedure."))
  if(redo){
    if(length(outlier.res$outidx) > 0){
      message("Redo mean estimation after removing outliers...\n")
      mean.mat = matrix(NA, M, N)
      if(is.null(cluster)){
        if(length(idx0) > 0)
          mean.mat[, !msk] = .pixelwiseMeanEst(doy[-idx0], mat[-idx0, !msk], doyeval, minimum.num.obs, NULL) else
            mean.mat[, !msk] = .pixelwiseMeanEst(doy, mat[, !msk], doyeval, minimum.num.obs, NULL)
      } else {
        if(length(idx0) > 0)
          mean.mat[, !msk] = .pixelwiseMeanEst(doy[-idx0], mat[-idx0, !msk], doyeval, minimum.num.obs, cluster[!msk]) else
            mean.mat[, !msk] = .pixelwiseMeanEst(doy, mat[, !msk], doyeval, minimum.num.obs, cluster[!msk])
      }
    }
  }
  
  ###########################################
  ######### Mean value correction ###########
  ###########################################
  ## Correct values that beyond clipRange or missing values in mean.mat
  message("Doing mean value correction...\n")
  if(clipMethod == "truncate"){
    mean.mat[mean.mat < clipRange[1]] = clipRange[1]
    mean.mat[mean.mat > clipRange[2]] = clipRange[2]
  } else
    if(clipMethod == "nnr"){
      mean.mat[mean.mat < clipRange[1] | mean.mat > clipRange[2]] = NA
      ## find columns that have NA values
      msk.idx = which(msk == 1)
      notmsk.idx = which(msk == 0)
      napixel.idx = which(apply(mean.mat, 2, function(x) any(is.na(x))))
      napixel.in.msk = napixel.idx %in% msk.idx
      napixel.idx = napixel.idx[!napixel.in.msk]
      if(length(napixel.idx) > 0){
        for(i in napixel.idx){
          d = 3
          ## neighbor pixel indexes
          if(is.null(img.nrow) | is.null(img.ncol))
            stop("img.nrow and img.ncol is not supplied, which are needed for clipMethod == nnr.")
          nbrpixel.idx = intersect(nbr(i-1, img.nrow, img.ncol, d, d) + 1, notmsk.idx)
          while(all(nbrpixel.idx %in% napixel.idx)){
            d = d + 2
            nbrpixel.idx = intersect(nbr(i-1, img.nrow, img.ncol, d, d) + 1, notmsk.idx)
          }
          ## mean of neighborhood pixels
          mm = apply(mean.mat[, nbrpixel.idx], 1, mean, na.rm=TRUE)
          mm.miss.idx = is.na(mean.mat[,i])
          mean.mat[mm.miss.idx, i] = mm[mm.miss.idx]
        }
      }
    }
  return(list(
    doyeval = doyeval,
    meanmat = mean.mat,
    idx = list(idx.allmissing = idx1,
               idx.partialmissing = idx2,
               idx.fullyobserved = idx3,
               idx.outlier = idx0),
    outlier = outlier.res
  ))
}

.pixelwiseMeanEst <- function(doy, mat, doyeval, minimum.num.obs, cluster){
  idx = 1:length(doy) ## idx is the index of image, 1, 2, 3,...
  temporal_mean_est = stfit::opts_stfit$get("temporal_mean_est")
  N = ncol(mat) ## number of pixels
  M = length(doyeval) ## number of doy to evaluate on
  if(length(doy) != nrow(mat))
    stop("doy and mat dimension do not match!")
  mean.mat = matrix(NA, M, N)
  if(is.null(cluster)){
    ## Estimate the overall mean curves for each pixel
    message("Estimating the overall mean curve for each pixel...\n")
    i = NULL
    mean.mat = foreach(i = 1:N, .combine = "cbind") %dopar% {
      temporal_mean_est(doy, mat[, i], doyeval, minimum.num.obs)
    }
  } else {
    if(!is.vector(cluster))
      stop("cluster should be a vector.")
    if(length(cluster) != N)
      stop("cluster dimension is not correct.")
    uc = unique(cluster) ## unique clusters w/o msk cluster
    ## Estimate the mean curves for each cluster
    message("Estimating mean curves for each cluster...\n")
    for(cl in uc){
      clidx = cluster == cl
      mean.mat[, clidx] = temporal_mean_est(rep(doy, sum(clidx)), c(mat[, clidx]), doyeval)
    }
  }
  return(mean.mat)
}

#' Smoothing spline regression
#'
#' @param x independent variable
#' @param y response variable
#' @param x.eval vector to predict on
#' @param minimum.num.obs minimum number of observations needed to run the regression
#' @param ... other parameters to be passed to \code{smooth.spline} function
#'
#' @return predicted values at 'x.eval'
#' @export
#' @example 
#' plot(cars$dist,cars$speed)
#' lines(1:120, smooth_spline(cars$dist, cars$speed, 1:120))
smooth_spline <- function(x, y, x.eval=x, minimum.num.obs=4, ...) {
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) > minimum.num.obs){
    splfit <- smooth.spline(x[nonna.idx], y[nonna.idx], ...)
    res = predict(splfit, x.eval)$y
    # res[x.eval < min(x) | x.eval > max(x)] = NA
    return(res)
  } else{
    return(rep(NA, length(x.eval)))
  }
}

#' Local linear regression
#'
#' @param x independent variable
#' @param y response variable
#' @param h bandwidth
#' @param Kern Kernel
#' @param x.eval dnew data to predict on
#' @param minimum.num.obs minimum number of observations needed to run the regression
#'
#' @return predicted values at 'x.eval'
#' @export
#' @example 
#' plot(cars$dist,cars$speed)
#' lines(1:120, llreg(cars$dist, cars$speed, 1:120))
llreg <- function(x, y, x.eval=x, minimum.num.obs = 4, h=60, Kern=epan){
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) > minimum.num.obs){
    res <- .llreg(x[nonna.idx], y[nonna.idx], x.eval, h, Kern)
    return(res)
  } else{
    return(rep(NA, length(x.eval)))
  }
}

.llreg <- function(x, y, x.eval, h, Kern){
  if(length(x) != length(y))
    stop("llreg: x and y should have the same length")
  N = length(x.eval)
  fitted = rep(0, N)
  for(i in 1:N){
    vecX = x - x.eval[i]
    vecK = Kern((vecX)/h)
    sumK = sum(vecK)
    h1 = h
    while(sumK == 0){
      h1 = h1 + 1
      vecK = Kern((vecX)/h1)
      sumK = sum(vecK)
    }
    S0 = mean(vecK)
    S1 = mean(vecX*vecK) ## higher order term than S0 and S2, keep this term for accuracy 
    S2 = mean(vecX^2*vecK)
    denom = S2*S0 - S1^2
    if(denom != 0){
      fitted[i] = mean((S2-S1*vecX)*vecK*y)/denom ## local linear estimator
    } else
      fitted[i] = sum(vecK*y)/sumK ## local constant estimator
  }
  fitted
}

#' Local Polynomial Regression
#'
#' @param x independent variable
#' @param y response variable
#' @param x.eval vector to predict on
#' @param minimum.num.obs minimum number of observations needed to run the regression
#' @param span see 'loess' function
#' @param ... other parameters passed to 'loess' function
#'
#' @return predicted values at 'x.eval'
#' @export
#' @example 
#' plot(cars$dist,cars$speed)
#' lines(1:120, lpreg(cars$dist, cars$speed, 1:120, span = 20))
lpreg <- function(x, y, x.eval, minimum.num.obs = 4, span=0.3, ...){
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) > minimum.num.obs){
    x = x[nonna.idx]
    y = y[nonna.idx]
    loessfit <- loess(y~x, span = span, control = loess.control(surface = "direct"), ...)
    res = predict(loessfit, data.frame(x = x.eval))
    return(res)
  } else{
    return(rep(NA, length(x.eval)))
  }
}

#' spline regression
#'
#' @param x independent variable
#' @param y response variable
#' @param x.eval vector to predict on
#' @param minimum.num.obs minimum number of observations needed to run the regression
#' @param basis what basis to use, "fourier" and "bspline" are available
#' @param rangeval see \code{fda::create.basis}
#' @param nbasis see \code{fda::create.basis}
#' @param ... arguments passed to \code{fad::create.basis} functions
#' @return predicted values at 'x.eval'
#' @export
#' @example 
#' plot(cars$dist,cars$speed)
#' lines(1:120, spreg(cars$dist, cars$speed, 1:120, basis = 'bspline', nbasis = 5))
spreg <- function(x, y, x.eval, minimum.num.obs = 4, basis = c("fourier", "bspline"), 
                  rangeval = c(min(x.eval)-1, max(x.eval)), nbasis = 11, ...){
  nonna.idx = !is.na(y)
  basis = match.arg(basis)
  if(sum(nonna.idx) > minimum.num.obs){
    x = x[nonna.idx]
    y = y[nonna.idx]
    # whisker = boxplot(y, plot = FALSE)$stats[c(1, 5)]
    # idx = (y > whisker[1]) & (y <whisker[2])
    # x = x[idx]
    # y = y[idx]
    if(basis == "fourier"){
      bs = fda::create.fourier.basis(rangeval=rangeval, nbasis=nbasis, ...)
    } else
      if(basis == "bspline"){
        bs = fda::create.bspline.basis(rangeval=rangeval, nbasis=nbasis, ...)
      }
    X = fda::eval.basis(x, bs)
    lmfit = lm.fit(X, y)
    return(fda::eval.basis(x.eval, bs) %*% lmfit$coefficients)
  } else{
    return(rep(NA, length(x.eval)))
  }
}


