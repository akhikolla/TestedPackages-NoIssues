#' This functions calculates numerical gradient of log-likelihood of cDCC-GARCH model.
#' @param param cDCC-GARCH parameters(alpha,beta)
#' @param ht matrix of conditional variance vectors (T by N)
#' @param residuals matrix of residual(de-mean) returns  (T by N)
#' @param stdresids matrix of standrdized(De-GARCH) residual returns  (T by N)
#' @param uncR unconditional correlation matrix of stdresids  (N by N)
#' @param d (log-lik(x+d) - log-lik(x))/d
#'
#' @return numerical gradient of log-likelihood of cDCC-GARCH model(vector)
#' @export

cdcc_gradient <- function(param, ht, residuals, stdresids, uncR, d=1e-5){
  nobs <- dim(residuals)[1]
  ndim <- dim(residuals)[2]
  npara <- length(param)
  Id <- d*diag(npara)
  param1 <- param + Id[,1]
  param2 <- param + Id[,2]

  lfc <- cdcc_compositelik(param[1], param[2], ht, residuals, stdresids,  uncR, nobs, ndim)
  lfc1 <- cdcc_compositelik(param1[1], param1[2], ht, residuals, stdresids,  uncR, nobs, ndim)
  lfc2 <- cdcc_compositelik(param2[1], param2[2], ht, residuals, stdresids,  uncR, nobs, ndim)

  c(sum((lfc1 - lfc)/d), sum((lfc2 - lfc)/d) )

}
