#' This function calculates log-likelihood of DCC-GARCH model.
#' @param param DCC-GARCH parameters(alpha,beta)
#' @param ht matrix of conditional variance vectors (T by N)
#' @param residuals matrix of residual(de-mean) returns  (T by N)
#' @param stdresids matrix of standrdized(De-GARCH) residual returns  (T by N)
#' @param uncR unconditional correlation matrix of stdresids  (N by N)
#'
#' @return log-likelihood of DCC-GARCH model (scaler)

dcc_loglikelihood <- function(param, ht, residuals, stdresids, uncR){
  nobs <- dim(stdresids)[1]
  ndim <- dim(stdresids)[2]

  lf <- dcc_compositelik(param[1], param[2], ht, residuals, stdresids,  uncR, nobs, ndim)

  return(lf)
}
