## Copyright 2011 – 2020  Center for Survey Statistics and Methodology at Iowa State University   All Rights Reserved
#' STFIT for Landsat data
#'
#' This function is used for Landsat data imputation, which includes five steps: 
#' mean estimation, outlier detection, temporal effect estimation, spatial effect
#' estimation and imputation. In real application, one can use this as a template
#' to create a five steps imputation procedure depending on the real data structure.
#'
#' @param year vecotr of year
#' @param doy vecotr of DOY (day of the year)
#' @param mat a numeric matrix. Each row contains a row stacked image pixel values.
#' @param img.nrow number of rows of the image
#' @param img.ncol number of columns of the image
#' @param doyeval a vector of DOY on which to get the mean and temporal imputation
#' @param h.tcov bandwidth for temporal covariance estimation
#' @param h.tsigma2 bandwith for temporal variance estimation 
#' @param h.scov bandwidth for spatial covariance estimation
#' @param h.ssigma2 bandwidth for spatial variance estimation
#' @param nnr maximum number of nearest neighbor pixels to use for spatial covariance estimation
#' @param outlier.action "keep" to keep outliers; "remove" to replace outliers with imputed values
#' @param intermediate.save TRUE or FASLE; whether to save the intermediate results including
#' mean, temporal effect and spacial effect imputation resutls. The intermediate results can be
#' useful to avoid duplicating the computation for some imputation steps.
#' @param intermediate.dir directory where to save the intermediate results
#' @param use.intermediate.result whether to use the intermediate results in the 'intermediate.dir' folder. 
#' Default is TRUE.
#' @param teff TRUE or FALSE, wheter to calculate the temporal effect. Default is TRUE.
#' @param seff TRUE or FALSE, wheter to calculate the spatial effect. Default is TRUE.
#' @param doy.break a vector of break points for \code{doy} where the spatial effect are 
#' estimated seperately on each interval. Default is NULL, i.e. the spatial effect is assumed
#' to be the same over \code{doy}.
#' @param cycle TRUE or FALSE. When \code{doy.break} is specified, whether to combine the first
#' \code{doy.break} interval and the last \code{doy.break} together for spatial effect estimation.
#' @param outlier.tol The threshold to use to define outlier image. Default is 0.2, i.e. images
#' with more than 20\% outlier pixels are treated as outlier image.
#' @param t.grid a vector of grid points on which to calculate the temporal covariance function
#' @param t.grid.num number of grid points to use for temporal covariance estimation. 
#' Ignored if \code{t.grid} is given.
#' @param clipRange passed to \code{meanEst} function
#' @param clipMethod passed to \code{meanEst} function
#' @param var.est Whether to estimate the variance of the temporal and spatial effects. Default is FALSE.
#'
#' @return List of length 4 with entries:
#' \itemize{
#'   \item imat: imputed matrix of \code{mat}
#'   \item smat: standard error matrix of the same size as \code{mat}
#'   \item idx: a list of image indexes 
#'   \itemize{
#'     \item idx.allmissing: completely missing image indexes,
#'     \item idx.partialmissing: partially observed image indexes,
#'     \item idx.fullyobserved: fully observed image indexes,
#'     \item idx.outlier: outlier image indexes.
#'   }
#'   \item outlier: a list of image outliers information
#'   \itemize{
#'     \item outidx: image index with outlier pixels,
#'     \item outpct: percentage of outlier pixels corresponding to \code{outidx},
#'     \item outlst: a list of the same length as \code{outidx}, with each list the missing 
#'     pixel index.
#'   }
#' }
#' @export
#'
#' @examples
#' \donttest{
#' library(doParallel)
#' library(raster)
#' library(rasterVis)
#' library(RColorBrewer)
#' dfB = landsat106[landsat106$year >= 2000,]
#' matB = as.matrix(dfB[,-c(1:2)])
#' year = dfB$year
#' doy = dfB$doy
#' if(require(doParallel))
#'   registerDoParallel(1)
#' res <- stfit_landsat(year, doy, matB, 31, 31, nnr=30,
#' use.intermediate.result = FALSE, intermediate.save = FALSE, var.est = TRUE)
#' ## visualize the imputed results
#' idx = c(res$idx$idx.allmissing[150], res$idx$idx.partialmissing[c(30, 60, 90)])
#' rst_list = list()
#' for(i in 1:length(idx)){
#'   rst_list[(i-1)*3+1] = raster(matrix(matB[idx[i],], 31))
#'   rst_list[(i-1)*3+2] = raster(matrix(res$imat[idx[i],], 31))
#'   rst_list[(i-1)*3+3] = raster(matrix(res$sdmat[idx[i],], 31))
#' }
#' s = stack(rst_list)
#' levelplot(s, index.cond=list(c(seq(1, 12, 3), seq(2, 12, 3), seq(3, 12, 3))),
#'           par.setting = rasterTheme(panel.background=list(col="black"),
#'                                     region = brewer.pal(9, 'YlOrRd')),
#'           names.attr = c(rbind(paste0("Original ", idx), 
#'                                paste0("Imputed ", idx),
#'                                paste0("Std. Error ", idx))),
#'           layout = c(4,3))
#' 
#' }
stfit_landsat <- function(year, doy, mat, img.nrow, img.ncol, doyeval = 1:365,  h.tcov = 100, h.tsigma2 = 300,
                            h.scov = 2, h.ssigma2 = 2, nnr = 10, outlier.action = c("keep", "remove"), outlier.tol = 0.2,
                            intermediate.save = TRUE, intermediate.dir = "./intermediate_output/",
                            use.intermediate.result = TRUE, teff = TRUE, seff = TRUE,
                            doy.break = NULL, cycle = FALSE, t.grid =NULL, t.grid.num = 50,
                            clipRange = c(0, 1800), clipMethod = "nnr", var.est = FALSE){
  if(intermediate.save){
    if(!dir.exists(intermediate.dir)){
      message(paste0("Folder", intermediate.dir, "is created to save intermediate results."))
      dir.create(intermediate.dir, recursive = TRUE)
    }
  }
  
  imat = mat
  ###################################
  #### 1. Overall mean estimaton ####
  ###################################
  if(use.intermediate.result & file.exists(paste0(intermediate.dir, "meanest.rds"))){
    meanest = readRDS(paste0(intermediate.dir, "meanest.rds"))
  } else {
    meanest = meanEst(doy, mat, doyeval = doyeval, outlier.tol = outlier.tol, clipRange = clipRange,
                      clipMethod = clipMethod, img.nrow = img.nrow, img.ncol = img.ncol)
    if(intermediate.save)
      saveRDS(meanest, paste0(intermediate.dir, "meanest.rds"))
  }
  
  ###################################
  #### 2. Time effect estimation ####
  ###################################
  ## remove outlier pixels
  for(i in 1:length(meanest$outlier$outidx)){
    mat[meanest$outlier$outidx[i], meanest$outlier$outlst[[i]]] = NA
  }
  ## remove outlier images
  outlier.img.idx = meanest$idx$idx.outlier
  for(i in outlier.img.idx){
    mat[outlier.img.idx,] = NA
  }
  ## calculate the residuals
  rmat = mat - meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]
  
  ## estimate the temporal effect using residuals
  ## result is a 3d array with the first dimension year, second dimension doy and third dimension pixel index
  if(teff){
    if(use.intermediate.result & file.exists(paste0(intermediate.dir, "teffres.rds"))){
      teffres = readRDS(paste0(intermediate.dir, "teffres.rds"))
    } else {
      teffres = teffEst(year, doy, rmat, doyeval = meanest$doyeval, h.cov = h.tcov, h.sigma2 = h.tsigma2, 
                        t.grid = t.grid, t.grid.num = t.grid.num, var.est = var.est)
      if(intermediate.save)
        saveRDS(teffres, paste0(intermediate.dir, "teffres.rds"))
    }
    
    ######################################
    #### 3. Spatial effect estimation ####
    ######################################
    ## claculate residuals after removing temporal effect
    yearidx = unlist(lapply(year, function(x, y)
      which(y == x), y = as.numeric(dimnames(teffres$teff_array)[[1]])))
    doyidx = unlist(lapply(doy, function(x, y)
      which(y == x), y = as.numeric(dimnames(teffres$teff_array)[[2]])))
    for (i in 1:nrow(rmat)) {
      rmat[i, ] = rmat[i, ] - teffres$teff_array[yearidx[i], doyidx[i], ]
    }
  }
  
  ## estimate the spatial effect using residuals
  ## result is a 3d array with the first dimension year, second dimension doy and third dimension pixel index
  if(seff){
    if(use.intermediate.result & file.exists(paste0(intermediate.dir, "seffres.rds"))){
      seffres = readRDS(paste0(intermediate.dir, "seffres.rds"))
      seff_mat = seffres$seff_mat
      seff_var_mat = seffres$seff_var_mat
    } else {
      if(is.null(doy.break)){
        seffres = seffEst(rmat, img.nrow, img.ncol, nnr = nnr, h.cov = h.scov, h.sigma2 = h.ssigma2, var.est = var.est)
        seff_mat = seffres$seff_mat
        seff_var_mat = seffres$seff_var_mat
      } else {
        seff_mat = matrix(0, nrow(rmat), ncol(rmat))
        if(var.est) 
          seff_var_mat[tmpidx,] = matrix(0, nrow(rmat), ncol(rmat))
        if(cycle){
          tmpidx = (doy <= doy.break[1]) | (doy > doy.break[length(doy.break)])
          seffres = seffEst(rmat[tmpidx,], img.nrow, img.ncol, 
                                                  nnr = nnr, h.cov = h.scov, h.sigma2 = h.ssigma2)
          seff_mat[tmpidx,] = seffres$seff_mat
          if(var.est)
            seff_var_mat[tmpidx,] = seffres$seff_var_mat
          
          for(i in 2:length(doy.break)){
            tmpidx = (doy <= doy.break[i]) & (doy > doy.break[i-1])
            seffres = seffEst(rmat[tmpidx,], img.nrow, img.ncol, 
                                       nnr = nnr, h.cov = h.scov, h.sigma2 = h.ssigma2)
            seff_mat[tmpidx,] = seffres$seff_mat
            if(var.est)
              seff_var_mat[tmpidx,] = seffres$seff_var_mat
          }
        } else {
          tmpidx = (doy <= doy.break[1])
          seffres = seffEst(rmat[tmpidx,], img.nrow, img.ncol, 
                            nnr = nnr, h.cov = h.scov, h.sigma2 = h.ssigma2)
          seff_mat[tmpidx,] = seffres$seff_mat
          if(var.est)
            seff_var_mat[tmpidx,] = seffres$seff_var_mat
          for(i in 2:length(doy.break)){
            tmpidx = (doy <= doy.break[i]) & (doy > doy.break[i-1])
            seffres = seffEst(rmat[tmpidx,], img.nrow, img.ncol, 
                              nnr = nnr, h.cov = h.scov, h.sigma2 = h.ssigma2)
            seff_mat[tmpidx,] = seffres$seff_mat
            if(var.est)
              seff_var_mat[tmpidx,] = seffres$seff_var_mat
          }
          if(doy.break[i] < max(doy)){
            tmpidx = (doy > doy.break[i])
            seffres = seffEst(rmat[tmpidx,], img.nrow, img.ncol, 
                              nnr = nnr, h.cov = h.scov, h.sigma2 = h.ssigma2)
            seff_mat[tmpidx,] = seffres$seff_mat
            if(var.est)
              seff_var_mat[tmpidx,] = seffres$seff_var_mat
          }
        }
      }
      if(intermediate.save)
        saveRDS(list(seff_mat = seff_mat, seff_var_mat = seff_var_mat), paste0(intermediate.dir, "seffres.rds"))
    }
  }
  
  #######################
  #### 4. Imputation ####
  #######################
  ## partially missing images: mean + time effect + spatial effect 
  ## all missing images: mean + time effect
  ## first calculate the theoretically imputed mat.
  mat_imputed = meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]
  if(teff){
    for(i in 1:nrow(mat_imputed)){
      mat_imputed[i,] = mat_imputed[i,] + teffres$teff_array[yearidx[i], doyidx[i],]
    }
  }
  if(seff){
    mat_imputed = mat_imputed + seff_mat
  }
  
  if(var.est){
    sdmat = matrix(0, nrow(mat), ncol(mat))
    if(teff){
      for(i in 1:nrow(sdmat)){
        sdmat[i,] = sdmat[i,] + teffres$teff_var_array[yearidx[i], doyidx[i],]
      }
    }
    if(seff){
      sdmat = sdmat + seff_var_mat
    }
    sdmat = sqrt(sdmat)
  } else
    sdmat = NULL
  
  #### final imputation
  outlier.action = match.arg(outlier.action)
  if(outlier.action == "keep"){
    sdmat[!is.na(imat)] = NA
    imat[is.na(imat)] = mat_imputed[is.na(imat)]
  } else if (outlier.action == "remove") {
    sdmat[!is.na(mat)] = NA
    imat = mat
    imat[is.na(imat)] = mat_imputed[is.na(imat)]
  }
  return(list(imat = imat, sdmat = sdmat, idx = meanest$idx, outlier = meanest$outlier))
}


