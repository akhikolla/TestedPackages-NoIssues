#' Fisher's rule for combining several p-values 
#'
#' Compute the negative of the log-sum for a vector of p-values.
#'
#'@details According to Fisher's rule, if the p-values are correlated, then this does not follow a simple chi-square 
#' mixture under the null.
#'
#'@param pv the vector of p-values to be combined together
#'
#'@return the Fisher combination of the p-values. See Details.
#'
#'
#'@export

comb_pvals <- function(pv){
  
  fisher <- -sum(log(pv))
  
  return(fisher)
}