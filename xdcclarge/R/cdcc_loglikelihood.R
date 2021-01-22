#' This function calculates log-likelihood of cDCC-GARCH model.
#' @param param cDCC-GARCH parameters(alpha,beta)
#' @param ht matrix of conditional variance vectors (T by N)
#' @param residuals matrix of residual(de-mean) returns  (T by N)
#' @param stdresids matrix of standrdized(De-GARCH) residual returns  (T by N)
#' @param uncR unconditional correlation matrix of stdresids  (N by N)
#'
#' @return log-likelihood of cDCC-GARCH model (scaler)
#' @export

cdcc_loglikelihood <- function(param, ht, residuals, stdresids, uncR){
  nobs <- dim(stdresids)[1]
  ndim <- dim(stdresids)[2]

  lf <- cdcc_compositelik(param[1], param[2], ht, residuals, stdresids,  uncR, nobs, ndim)

  return(lf)
}
