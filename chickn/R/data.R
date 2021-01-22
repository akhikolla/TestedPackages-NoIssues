#' @title UPS2 dataset 
#' @description Proteomics data acquired within the mass spectrometry analysis of UPS2 sample.
#' @details Only a small part of data was taken from 
#' the original dataset described in \insertCite{henning2018peptide}{chickn}.  
#' The UPS2 dataset contains 190 chromatographic traces (matrix columns) acquired along 
#' the retention time (matrix rows) in the liquid chromatography.  
#' @source \url{https://github.com/optimusmoose/ups2GT} 
#' @references
#' \itemize{
#' \item \insertRef{tsou2015dia}{chickn} 
#' \item \insertRef{henning2018peptide}{chickn}
#' }
#' @format A matrix with 1653 rows and 190 columns. 
#' @usage data(UPS2)
"UPS2"