#' @title Positive Definite Quadratic Forms Distribution
#'
#' @description Density function, distribution function, quantile function and random generator for positive definite QFs.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param obj \code{MellinQF} object produced by the \code{\link{compute_MellinQF}} function.
#' @param eps_quant relative error for quantiles.
#' @param maxit_quant maximum number of Newton-Raphson iterations allowed to compute quantiles.
#' @param lambdas vector of positive weights.
#' @param etas vector of non-centrality parameters. Default all zeros.
#'
#'
#' @details
#'The quadratic form CDF and PDF are evaluated by numerical inversion of the Mellin transform.
#'The absolute error specified in \code{\link{compute_MellinQF}} is guaranteed for values of \code{q} and \code{x} inside the \code{range_q}.
#'If the quantile is outside \code{range_q}, computations are carried out, but a warning is sent.
#'
#'The function uses the Newton-Raphson algorithm to compute the QF quantiles related to probabilities \code{p}.
#'
#' @return
#' \code{dQF} provides the values of the density function at a quantile \code{x}.
#'
#' \code{pQF} provides the cumulative distribution function at a quantile \code{q}.
#'
#' \code{qQF} provides the quantile corresponding to a probability level \code{p}.
#'
#'  \code{rQF} provides a sample of \code{n} independent realizations from the QF.
#'
#'@seealso
#'See \code{\link{compute_MellinQF}} for details on the Mellin computation.
#'
#'
#' @examples
#'\donttest{
#' library(QF)
#' # Definition of the QF
#' lambdas_QF <- c(rep(7, 6),rep(3, 2))
#' etas_QF <- c(rep(6, 6), rep(2, 2))
#' # Computation Mellin transform
#' eps <- 1e-7
#' rho <- 0.999
#' Mellin <- compute_MellinQF(lambdas_QF, etas_QF, eps = eps, rho = rho)
#' xs <- seq(Mellin$range_q[1], Mellin$range_q[2], l = 100)
#' # PDF
#' ds <- dQF(xs, Mellin)
#' plot(xs, ds, type="l")
#' # CDF
#' ps <- pQF(xs, Mellin)
#' plot(xs, ps, type="l")
#' # Quantile
#' qs <- qQF(ps, Mellin)
#' plot(ps, qs, type="l")
#' #Comparison computed quantiles vs real quantiles
#' plot((qs - xs) / xs, type = "l")
#'
#'}
#'
#'@name QF

NULL




#' @name dQF
#' @rdname QF
#' @export

dQF <- function(x, obj) {
  if(!(class(obj)=="MellinQF" )){
    stop("A 'MellinQF' object is required for 'obj'")
  }
  if(sum((x < obj$range_q[1]) | (x > obj$range_q[2])) != 0) {
    warning("At least one element of 'x' is outside the range for which the desired error is assured")
  }
  .Call(`_QF_dQF_c_scal`, x, obj$Mellin)
}


#' @name pQF
#' @rdname QF
#' @export

pQF <- function(q, obj) {
  if(!(class(obj)=="MellinQF")){
    stop("A 'MellinQF' object is required for 'obj'")
  }
  if(sum((q < obj$range_q[1]) | (q > obj$range_q[2])) != 0) {
    warning("At least one element of 'q' is outside the range for which the desired error is assured")
  }
  .Call(`_QF_pQF_c_scal`, q, obj$Mellin)
}

#' @name qQF
#' @rdname QF
#' @export

qQF <- function(p, obj, eps_quant = 1e-6, maxit_quant = 1e4) {
  if(!(class(obj)=="MellinQF")){
    stop("A 'MellinQF' object is required for 'obj'")
  }
  if(sum((p < 0) | (p > 1)) != 0) {
    stop("At least one element of 'q' is outside the range (0,1)")
  }
  if(sum(eps_quant<0 || maxit_quant<0) != 0) {
    stop(" 'eps_quant' and 'maxit_quant' must be positive")
  }
  if(sum((p < (1 - obj$rho) / 2) | (p > 1-(1 - obj$rho) / 2)) != 0) {
    warning("At least one element of 'p' is outside the range for which the desired error is assured")
  }
  .Call(`_QF_qQF_c`, p, obj$Mellin, eps_quant,
        maxit_quant, sum(obj$lambdas)/obj$Mellin$lambda_min) * obj$Mellin$lambda_min

}


#' @name rQF
#' @rdname QF
#' @export



rQF <- function(n, lambdas, etas = rep(0, length(lambdas))) {
  .Call(`_QF_rQF_c`, n, lambdas, etas)
}

