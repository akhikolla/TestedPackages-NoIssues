#' Estimate Parameters of Dirichlet Distribution
#'
#' C++ implementation of the fixed-point iteration algorithm by Minka (2000).
#'
#' @param x a matrix of Dirichlet samples, one row per observation.
#' @param const constant that is added to avoid problems with zeros in \code{log(x)}.
#'     The default is \code{const = min(x[x>0])*.01}.
#' @param maxit maximum number of iterations.
#' @param abstol The absolute convergence tolerance: maximum of absolute differences of Dirichlet parameters.
#'
#' @details The algorithm is used to estimate the effective sample size based on samples
#'    of posterior model probabilities (see \code{\link{stationary}} and
#'    \code{\link{summary.stationary}}).
#'
#' @examples
#' x <- rdirichlet(100, c(8,1,3,9))
#' fit_dirichlet(x)
#'
#' @seealso \code{\link{rdirichlet}}
#' @references
#' Minka, T. (2000). Estimating a Dirichlet distribution. Technical Report.
#' @export
fit_dirichlet <- function (x, const, maxit = 1e5, abstol = .1){
  # adjust for x=0  (because of log(x) = -Inf)
  if (min(x) == 0){
    if (missing(const))
      const <- min(x[x>0])*.01
    x <- (x + const)/(1 + 2 * const)
  }
  x <- x/rowSums(x)
  logx.mean <- colMeans(log(x))
  N <- nrow(x)

  # heuristic for starting values:
  x.mean <- colMeans(x)
  x.squares <- colMeans(x^2)
  xi <- (x.mean - x.squares)/(x.squares - x.mean^2)
  alpha0 <- xi * x.mean
  alpha <- NULL
  try(alpha <- dirichlet_fp(pmax(.01, alpha0), logx.mean, # pmin(alpha0, 50)
                        maxit = maxit, abstol = abstol), silent = TRUE)
  # if this fails: random starting values
  if (is.null(alpha) || anyNA(alpha))
    alpha <- dirichlet_fp(runif(length(alpha),0.5,1),
                          logx.mean,
                          maxit = maxit, abstol = abstol)

  list("alpha" = alpha, "sum" = sum(alpha))
}
