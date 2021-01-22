## Copyright 2011 – 2020  Center for Survey Statistics and Methodology at Iowa State University   All Rights Reserved
#' stfit: Spatial-Temporal Functional Imputation Tool
#'
#' The stfit package provides functions to impute missing values for a 
#' sequence of observed images for the same location using functional
#' data analysis technique
#' @useDynLib stfit, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom foreach %dopar% foreach
#' @importFrom rasterVis rasterTheme levelplot 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom raster stack raster nlayers values calc
#' @importFrom stats lm.fit loess loess.control predict setNames smooth.spline var
#' @importFrom graphics boxplot
#' @importFrom Matrix sparseMatrix
#' @docType package
#' @name stfit-package
NULL
