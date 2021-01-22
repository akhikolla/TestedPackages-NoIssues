#' Generate a sample of a discrete-state Markov chain
#'
#' Generates a sequence of discrete states from a discrete-time Markov chain with transition matrix \code{P}.
#'
#' @param n length of the generated sequence.
#' @param P transition matrix (rows are normalized to sum to 1).
#' @param start vector with nonnegative values that defines the discrete starting
#'     distribution at t=0 (\code{start} is normalized to sum to 1). The default is
#'     a discrete uniform distribution.
#' @examples
#' P <- matrix(c(.30, .50, .20,
#'               .05, .25, .70,
#'               .00, .10, .90), 3, 3, byrow=TRUE)
#' rmarkov(50, P)
#' @export
rmarkov <- function(n, P, start = rep(1, ncol(P))){

  if(!is.matrix(P) || ncol(P) != nrow(P) ||
     any(P<0) || ncol(P) < 2)
      stop("'P' must be a transition matrix with positive values.")
  P <- P/rowSums(P)
  if(!is.numeric(n) || n != round(n) || n <= 1 || length(n) != 1)
    stop("'n' must be a positive integer")
  if(!is.numeric(start) || !is.vector(start) || length(start) != ncol(P))
    stop("Length of 'start' does not match size of 'P'.")

  c0 <- sample(1:ncol(P), 1, prob=start)
  z <- sim_mc(n, P, c0)
  c(z)
}
