#' MLE for the Zero-truncated Poisson Lognormal mixture distribtuion
#'
#' `ztplnmMLE` fits the Zero-truncated Poisson lognormal mixture distribution
#' to data and estimates parameters mean `mu`, standard deviation `sig` and 
#' mixture weight `theta` in the lognormal distribution.
#'
#' The function searches the maximum likelihood estimators of mean vector `mu`,
#' standard deviation vector `sig` and mixture weight vector `theta` using the
#' optimization procedures in \code{\link{nlminb}}.
#'
#' @param n a vector of counts
#' @param K number of components
#' @param lower_mu,upper_mu numeric values of lower and upper bounds for mean of
#' the variables's natural logarithm.
#' @param lower_sig,upper_sig numeric values of lower and upper bounds for
#' standard deviation of the variables's natural logarithm
#' @param lower_theta,upper_theta numeric values of lower and upper bounds for
#'  mixture weights.
#' @param type1 logical; if TRUE, Use type 1 ztpln else use type 2.  
#' @param message mean of lognormal distribution in sample 3.
#' @return \item{convergence}{An integer code. 0 indicates successful
#' convergence.}
#' @return \item{iterations}{Number of iterations performed.}
#' @return \item{message}{A character string giving any additional information
#' returned by the optimizer, or NULL. For details, see PORT documentation.}
#' @return \item{evaluation}{Number of objective function and gradient function
#'  evaluations}
#' @return \item{mu}{Maximum likelihood estimates of mu}
#' @return \item{sig}{Maximum likelihood estimates of sig}
#' @return \item{theta}{Maximum likelihood estimates of theta}
#' @return \item{loglik}{loglikelihood}
#' @examples
#' y <- rztplnm(100, c(1, 10), c(2, 1), c(0.2, 0.8))
#' ztplnmMLE(y)
#' @export
ztplnmMLE <- function(n,
                      K = 2,
                      lower_mu = rep(0, K),
                      upper_mu = rep(log(max(n)), K),
                      lower_sig = rep(0.001, K),
                      upper_sig = rep(10, K),
                      lower_theta = rep(0.001, K),
                      upper_theta = rep(0.999, K),
                      type1 = TRUE,
                      message = FALSE) {

  if (K == 1) {warning("`ztplnMLE` is faster for K = 1")}

  n <- n[n > 0]
  max_abund <- log(max(n))

  if (message) {
    sink("/dev/null")  # now suppresses
    fit0 <- try(mixtools::normalmixEM(log(n), k = K, maxrestart = 20))
    sink()
  } else {
    fit0 <- try(mixtools::normalmixEM(log(n), k = K, maxrestart = 20))
  }

  # set initial values using normalmixEM for better convergence
  if (class(fit0) == "try-error") {
    params <- c(rep(1, 2 * K), rep(1 / K, K))
  } else {
    params <- c(fit0$mu, fit0$sig, fit0$lambda)
  }

  loglik <- function(params, n) {
    mu <- params[1:K]
    sig <- params[(K + 1):(2 * K)]
    theta <- params[(2 * K + 1):(3 * K)]
    theta <- theta / sum(theta)
    if (type1) {
      lik <- do_dztplnm(n, mu, sig, theta)
    } else lik <- do_dztplnm2(n, mu, sig, theta)
    return(-sum(log(lik), na.rm = TRUE))
  }

  fit <- try(nlminb(start = params, # inital parameters
                    loglik, # function to be minimized (i.e., -logLike)
                    gradient = NULL,
                    lower = c(lower_mu, lower_sig, lower_theta),
                    upper = c(upper_mu, upper_sig, upper_theta),
                    n = n,
                    TRUE))
  
  fit$mu <- fit$par[1:K]
  fit$sig <- fit$par[(K + 1):(2 * K)]
  theta <- fit$par[(2 * K + 1):(3 * K)]
  fit$theta <- theta / sum(theta)
  fit$par <- NULL
  fit$loglik <- -fit$objective
  fit$objective <- NULL
  return(fit)
}
