## Copyright 2011 – 2020  Center for Survey Statistics and Methodology at Iowa State University   All Rights Reserved
#' Data visualization for landsat data
#'
#' @param mat A matrix, each row corresponds to a vectorized image pixel values.
#' @param img.nrow number of rows of the image
#' @param byrow logical value indicating whether the pixcel values are stored 
#' by row or by column. Default to FALSE
#' @param colthm Color theme for the plot, passing to the \code{par.settings} 
#' parameter of the \code{levelplot} function in the \code{rasterVis} package
#' @param ... All other options passed to  \code{levelplot} function in the 
#' \code{rasterVis} package
#'
#' @export
#'
#' @examples
#' landsatVis(landsat106[landsat106$year == 2015, -c(1:2)], 
#' names.attr = as.character(landsat106$doy[landsat106$year == 2015]))
landsatVis <- function(mat, img.nrow = 31, byrow=FALSE,
                       colthm = rasterTheme(panel.background=list(col="black"), 
                                            region = brewer.pal(9, 'YlOrRd')), ...){
  mat = as.matrix(mat)
  if(ncol(mat) %% img.nrow != 0)
    stop("Image dimension does not seem to be right: ncol(mat) %% img.nrow != 0")
  s = stack(lapply(1:nrow(mat), 
                   function(i) raster(matrix(mat[i,], nrow = img.nrow, byrow = byrow))))
  levelplot(s, par.settings = colthm, ...)
}

