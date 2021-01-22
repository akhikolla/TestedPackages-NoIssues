## Copyright 2011 – 2020  Center for Survey Statistics and Methodology at Iowa State University   All Rights Reserved
#' STFIT Spatial Effect Estimation
#' 
#' @param rmat residual matrix
#' @param img.nrow image row dimension
#' @param img.ncol image column dimension
#' @param pve percent of variance explained of the selected eigen values. Default is 0.99.
#' @param h.cov bandwidth for spatial covariance estimation; ignored if \code{weight.cov} is supplied
#' @param h.sigma2 bandwidth for sigma2 estimation
#' @param weight.cov weight matrix for spatial covariance estimation
#' @param weight.sigma2 weight vector for spatial variance estimation
#' @param nnr maximum number of nearest neighbor pixels to use for spatial covariance estimation
#' @param method "lc" for local constant covariance estimation and "emp" for empirical covariance estimation
#' @param msk an optional logistic vector. TRUE represent the corresponding pixel is always missing.
#' @param msk.tol if 'msk' is not given, the program will determine the mask using \code{getMask}
#' function. If the percentage of missing values for a pixel over time is greater than this
#' @param partial.only calculate the spatical effect for partially observed images only, default is TRUE
#' @param var.est Whether to estimate the variance of the temporal effect. Default is FALSE.
#' @return List of length 3 with entries:
#'   \itemize{
#'     \item seff_mat: estimated spatial effect matrix of the same shape as \code{rmat}.
#'     \item seff_var_mat: estimated spatial effect variance matrix of the same shape as \code{rmat}.
#'     \item idx: a list of two entries:
#'     \itemize{
#'       \item idx.allmissing: index of the completely missing images.
#'       \item idx.imputed: index of the partially observed images, where spatial effects 
#'       are estimated.
#'     }
#'   }
#' 
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
#' seffres = seffEst(rmat, 31, 31, nnr = 30)
#' idx = c(1:nrow(rmat))[apply(rmat, 1, function(x) {!all(is.na(x)) & sum(is.na(x)) != 0})]
#' landsatVis(seffres$seff_mat[idx[1:16],])
#' }
seffEst <- function(rmat, img.nrow, img.ncol, h.cov = 2, h.sigma2 = 2,
                    weight.cov = NULL, weight.sigma2 = NULL,
                    nnr, method = c("lc", "emp"), partial.only = TRUE,
                    pve = 0.99, msk = NULL, msk.tol = 0.95, var.est = FALSE){
  message("Estimating spatial effect...")
  if(is.null(weight.cov))
    weight.cov = weightMatrix(h.cov)
  if(is.null(weight.sigma2))
    weight.sigma2 = weightVector(h.sigma2)
  ## initialize the rmat imputed matrix
  idx = 1:nrow(rmat)
  N = img.nrow * img.ncol ## total number of pixels (including pixels in mask)
  if(ncol(rmat) != N)
    stop("number of columns of rmat does not match img.nrow*img.ncol.")

  if(is.null(msk)){
    msk = getMask(rmat, tol = msk.tol) # idx for 'black holes'
  } else {
    if(!is.vector(msk))
      stop("msk should be a vector.")
    if(length(msk) != N)
      stop("msk dimension is not correct.")
  }
  ## number of 'actual' pixels for each image (except for pixels in mask)
  N1 = sum(!msk) 
  
  ## initial spatial effect matrix
  seff_mat  = matrix(0, nrow(rmat), ncol(rmat))
  seff_mat[,msk] = NA
  seff_var_mat = seff_mat
  
  ## no spatial effect if there are too few actual pixels in one image
  if(N1 <= 1){
    warning('Too few pixels in the image')
    return(list(seff_mat=seff_mat,
                seff_var_mat=seff_var_mat,
                idx = list()))
  }
   
  pidx = (1:N)[!msk] ## 'actual' pixel indexes
  rmat = rmat[,!msk, drop=FALSE] ## now 'rmat' has N1 columns
  
  ## remove images that have 100% missing pixels
  pct_missing = apply(rmat, 1, function(x) {sum(is.na(x))/N1})
  ## Index
  idx1 = idx[pct_missing == 1] ## all missing images indexes;
  if(partial.only){
    idx2 = idx[pct_missing > 0 & pct_missing < 1] # partially missing images indexes;
  } else {
    idx2 = idx[pct_missing < 1] # partially missing + fully observed images indexes;
  }
  
  ###################################################
  ######### Covariance matrix estimation ############
  ###################################################
  message("Estimating the covariance matrix...\n")
  ## using fully observed + partially observed images for covariance estimation
  method <- match.arg(method)
  if(method == "lc"){
    if(N == N1)
      covest = sparse_lc_cov_est(rmat, weight.cov, img.nrow, img.ncol, nnr)
    else
      covest = sparse_lc_cov_est1(rmat, weight.cov, img.nrow, img.ncol, nnr, pidx - 1)
  } else
    if (method == "emp") {
      if (N == N1)
        ## no black hole
        covest = sparse_emp_cov_est(rmat, img.nrow, img.ncol, nnr)
      else
        covest = sparse_emp_cov_est1(rmat, img.nrow, img.ncol, nnr, pidx -
                                       1) #pidx start from 0 in C++
    }
  ## coerce NA to 0; two vectors with no overlapped observations leads to NA
  covest$value[is.na(covest$value)] = 0
  scovest = sparseMatrix(covest$ridx, covest$cidx, x = covest$value, 
                         dims = c(N1, N1), symmetric = TRUE)
  message("Estimating the variance function...\n")
  ## using empirical variance estimation; TODO: using local constant/linear smoothing
  sigma2 = apply(rmat, 2, var, na.rm=TRUE)
  nugg = max(0, mean(sigma2 - Matrix::diag(scovest), na.rm=TRUE))
  
  message("Doing eigen decomposition on covariance matrix...\n")
  ## Eigen decomposition of the covariance matrix; may be improved later using rARPACK
  ev = eigen(scovest)
  ev$values[ev$values < 0] = 0
  
  ev.idx = which.min(cumsum(ev$values)/sum(ev$values) < pve)
  message("Spatial effect covariance estimation: the first ", ev.idx, " eigen values are used...\n")
  ev.vec = ev$vectors[, 1:ev.idx, drop=FALSE]
  ev.val = ev$values[1:ev.idx]
  
  ###############################################
  ######### Missing value imputation ############
  ###############################################
  ## The following uses the sparse pca to impute the missing values.
  message("Estimating the principal component scores for partially missing images...\n")
  seffres = PACE2d(rmat[idx2,, drop = FALSE], ev.vec, nugg, ev.val, var.est = var.est)
  seff_mat[idx2, !msk] = seffres$smat
  if(var.est)
    seff_var_mat[idx2, !msk] = seffres$svarmat 
  else
    seff_var_mat = NULL

  return(list(
    seff_mat = seff_mat,
    seff_var_mat = seff_var_mat,
    idx = list(idx.allmissing = idx1,
               idx.imputed = idx2)
  ))
  
}




