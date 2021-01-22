#' MLE for the Zero-truncated Poisson Lognormal distribution
#'
#' `ztplnMLE` fits the Zero-truncated Poisson lognormal distribution to data and
#'  estimates parameters mean `mu` and standard deviation `sig` in the lognormal
#'  distribution
#'
#'  The function searches the maximum likelihood estimates of mean `mu` and
#'  standard deviation `sig` using the optimization procedures in
#'  \code{\link{nlminb}}.
#'
#' @param n a integer vector of counts
#' @param lower_mu,upper_mu numeric values of lower and upper bounds for mean of
#'  the variables's natrual logarithm.
#' @param lower_sig,upper_sig numeric values of lower and upper bounds for
#' standard deviatoin of the variables's natrual logarithm
#' @param type1 logical; if TRUE, Use type 1 ztpln else use type 2.  
#' @return \item{convergence}{An integer code. 0 indicates successful
#' convergence.}
#' @return \item{iterations}{Number of iterations performed.}
#' @return \item{message}{A character string giving any additional information
#' returned by the optimizer, or NULL. For details, see PORT documentation.}
#' @return \item{evaluation}{Number of objective function and gradient function
#'  evaluations}
#' @return \item{mu}{Maximum likelihood estimates of mu}
#' @return \item{sig}{Maximum likelihood estimates of sig}
#' @return \item{loglik}{loglikelihood}
#' @examples
#' y <- rztpln(100, 3, 2)
#' ztplnMLE(y)
#' @export
ztplnMLE <- function(n,
                     lower_mu = 0,
                     upper_mu = log(max(n)),
                     lower_sig = 0.001,
                     upper_sig = 10,
                     type1 = TRUE) {

  n <- n[n > 0]
  # set initial values using medians and sd for fast convergence
  params <- c(log(median(n)), sd(log(n)))

  loglik <- function(params, n) {
    mu <- params[1]
    sig <- params[2]
    if (type1) {
      lik <- do_dztpln(n, mu, sig)
    } else lik <- do_dztpln2(n, mu, sig)
    return(-sum(log(lik), na.rm = TRUE))
  }
  # faster than Nelder-Mead and BFGS in optim
  fit <- try(nlminb(params, # initial parameters
                    loglik, # function to be minimized (i.e., -logLike)
                    gradient = NULL,
                    lower = c(lower_mu, lower_sig),
                    upper = c(upper_mu, upper_sig),
                    n = n,
                    TRUE))
  fit$mu <- fit$par[1]
  fit$sig <- fit$par[2]
  fit$par <- NULL
  fit$loglik <- -fit$objective
  fit$objective <- NULL
  return(fit)
}
