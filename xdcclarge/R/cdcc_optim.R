#' This function optimizes log-likelihood of cDCC-GARCH model.
#' @param param cDCC-GARCH parameters(alpha,beta)
#' @param ht matrix of conditional variance vectors (T by N)
#' @param residuals matrix of residual(de-mean) returns  (T by N)
#' @param stdresids matrix of standrdized(De-GARCH) residual returns  (T by N)
#' @param uncR unconditional correlation matrix of stdresids  (N by N)
#'
#' @return results of optimization
#' @export

cdcc_optim <- function(param, ht, residuals, stdresids, uncR){

  resta <- rbind(c(-1, -1), diag(2))
  restb <- c(-1, 0, 0)

  opt <- stats::constrOptim(theta=param, f=cdcc_loglikelihood, grad=cdcc_gradient, ui=resta, ci=restb, mu=1e-5, outer.iterations=10, outer.eps=1e-2, control=list(maxit=10,reltol=1e-5), ht=ht, residuals=residuals, stdresids=stdresids, uncR=uncR)

  return(opt)

}
