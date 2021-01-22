#' Compute p-values for a Z-score
#' 
#' Compute p-values for a Z-score assuming normal distribution of the z-score 
#' under the null Hypothesis H0
#' 
#' @param beta the estimate
#' @param sigma estimate's estimated variance
#' 
#' @return the p-value
#' 
#' @importFrom stats pnorm
#' 
#' @export

pval_zscore <- function(beta, sigma){
  
  2*stats::pnorm(-abs(beta/sigma))
  
}