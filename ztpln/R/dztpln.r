#' The zero-truncated compund poisson-lognormal distributions
#'
#' Density function and random generation for Zero-Trauncated Poisson Lognormal
#' distribution with parameters `mu` and sd `sig`.
#'
#' A compound Poisson-lognormal distribution is a Poisson probability 
#' distribution where its parameter \eqn{\lambda} is a random variable with 
#' lognormal distribution, that is to say \eqn{log\lambda} are normally 
#' distributed with mean \eqn{\mu} and variance \eqn{\sigma^2} (Bulmer 1974). 
#' The zero-truncated Poisson-lognormal distribution can be derived from a 
#' zero-truncated Poisson distribution.
#' 
#' Type 1 ZTPLN truncates zero based on Poisson-lognormal distribution and 
#' type 2 ZTPLN truncates zero based on zero-truncated Poisson distribution.
#' For mathematical details, please see `vignette("ztpln")`
#'
#' @param n number of random values to return.
#' @param x	vector of (non-negative integer) quantiles.
#' @param mu mean of lognormal distribution.
#' @param sig standard deviation of lognormal distribution.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param type1 logical; if TRUE, Use type 1 ztpln else use type 2.  
#' @return dztpln gives the (log) density and rztpln generates
#' random variates.
#' @references Bulmer, M. G. 1974. On Fitting the Poisson Lognormal Distribution to Species-Abundance Data. Biometrics 30:101-110.
#'
#' @seealso \code{\link{dztplnm}}
#'
#' @examples
#' rztpln(n = 10, mu = 0, sig = 1, type1 = TRUE)
#' rztpln(n = 10, mu = 6, sig = 4, type1 = TRUE)
#' dztpln(x = 1:5, mu = 1, sig = 2)
#' @export
dztpln <- function(x, mu, sig, log = FALSE, type1 = TRUE) {
  if (length(mu) > 1 | length(sig) > 1)
    stop("Vectorization of parameters not implemented")
  if (sig < 0) stop("sig needs to be > 0")
  if (any(!DistributionUtils::is.wholenumber(x))) warning("non integer values in x")
  if (min(x) <= 0) warning("zero in x")
  x <- x[x > 0]
  if (type1) {
    lik <- do_dztpln(x, mu, sig)
  } else lik <- do_dztpln2(x, mu, sig)
  if (log) return(log(lik)) else return(lik)
}

#' @rdname dztpln
#' @export
rztpln <- function(n, mu, sig, type1 = TRUE) {
  if (length(mu) > 1 | length(sig) > 1) 
     stop("Vectorization of parameters not implemented")
  if (sig < 0) stop("sig needs to be > 0")
  if (sig > 15) stop("standard deviation > 15 in log-scale is too large")
  if (type1) do_vec_rztpln1(n, mu, sig) else do_vec_rztpln2(n, mu, sig)
}
