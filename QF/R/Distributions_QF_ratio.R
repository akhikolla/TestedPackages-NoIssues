#' @title Ratio of Positive Definite Quadratic Forms Distribution
#'
#' @description Density function, distribution function, quantile function and random generator for the ratio of positive definite QFs.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param obj \code{MellinQF_ratio} object produced by the \code{\link{compute_MellinQF_ratio}} function.
#' @param eps_quant relative error for quantiles.
#' @param maxit_quant maximum number of Newton-Raphson iterations allowed to compute quantiles.
#' @param lambdas_num vector of positive weights of the QF at the numerator.
#' @param lambdas_den vector of positive weights of the QF at the denominator
#' @param etas_num vector of non-centrality parameters of the QF at the numerator. Default all zeros.
#' @param etas_den vector of non-centrality parameters of the QF at the denominator Default all zeros.
#'
#'
#' @details
#'The CDF and PDF of the ratio of positive QFs are evaluated by numerical inversion of the Mellin transform.
#'The absolute error specified in \code{\link{compute_MellinQF_ratio}} is guaranteed for values of \code{q} and \code{x} inside \code{range_q}.
#'If the quantile is outside \code{range_q}, computations are carried out, but a warning is sent.
#'
#'The function uses the Newton-Raphson algorithm to compute the ratio of QFs quantiles related to probabilities \code{p}.
#'
#' @return
#' \code{dQF_ratio} provides the values of the density function at a quantile \code{x}.
#'
#'\code{pQF_ratio} provides the cumulative distribution function at a quantile \code{q}.
#'
#' \code{qQF_ratio} provides the quantile corresponding to a probability level \code{p}.
#'
#' \code{rQF_ratio} provides a sample of \code{n} independent realizations the QFs ratio.
#'
#'@seealso
#'See \code{\link{compute_MellinQF_ratio}} for details on the Mellin computation.
#'
#'
#' @examples
#'\donttest{
#' lambdas_QF_num <- c(rep(7, 6),rep(3, 2))
#' etas_QF_num <- c(rep(6, 6), rep(2, 2))
#' lambdas_QF_den <- c(0.6, 0.3, 0.1)
#' # Computation Mellin transform
#' eps <- 1e-7
#' rho <- 0.999
#' Mellin_ratio <- compute_MellinQF_ratio(lambdas_QF_num, lambdas_QF_den,
#'                                        etas_QF_num, eps = eps, rho = rho)
#' xs <- seq(Mellin_ratio$range_q[1], Mellin_ratio$range_q[2], l = 100)
#' # PDF
#' ds <- dQF_ratio(xs, Mellin_ratio)
#' plot(xs, ds, type="l")
#' # CDF
#' ps <- pQF_ratio(xs, Mellin_ratio)
#' plot(xs, ps, type="l")
#' # Quantile
#' qs <- qQF_ratio(ps, Mellin_ratio)
#' plot(ps, qs, type="l")
#' #Comparison computed quantiles vs real quantiles
#' plot((qs - xs) / xs, type = "l")
#'
#'}
#'
#'@name QF_ratio

NULL


#' @name dQF_ratio
#' @rdname QF_ratio
#' @export

dQF_ratio <- function(x, obj) {
  if(! class(obj)=="MellinQF_ratio"){
    stop("A 'MellinQF_ratio' object is required for 'obj'")
  }
  if(sum((x < obj$range_q[1]) | (x > obj$range_q[2])) != 0) {
    warning("At least one element of 'x' is outside the range for which the desired error is assured")
  }
  .Call(`_QF_dQF_c_scal`, x, obj$Mellin)
}


#' @name pQF_ratio
#' @rdname QF_ratio
#' @export

pQF_ratio <- function(q, obj) {
  if(! class(obj) == "MellinQF_ratio"){
    stop("A 'MellinQF_ratio' object is required for 'obj'")
  }
  if(sum((q < obj$range_q[1]) | (q > obj$range_q[2])) != 0) {
    warning("At least one element of 'q' is outside the range for which the desired error is assured")
  }
  .Call(`_QF_pQF_c_scal`, q, obj$Mellin)
}

#' @name qQF_ratio
#' @rdname QF_ratio
#' @export

qQF_ratio <- function(p, obj, eps_quant = 1e-6, maxit_quant = 1e4) {
  if(! class(obj) == "MellinQF_ratio"){
    stop("A 'MellinQF_ratio' object is required for 'obj'")
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
        maxit_quant, (obj$range_q[1]+obj$range_q[2])/(2*obj$Mellin$lambda_min)) * obj$Mellin$lambda_min
}





#' @name rQF_ratio
#' @rdname QF_ratio
#' @export


rQF_ratio <- function(n, lambdas_num, lambdas_den,
                      etas_num = rep(0, length(lambdas_num)),
                      etas_den = rep(0, length(lambdas_den))) {
  .Call(`_QF_rQF_ratio_c`, n, lambdas_num, lambdas_den, etas_num, etas_den)
}


