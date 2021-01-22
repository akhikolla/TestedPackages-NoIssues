#' @title
#' Probability weighted moments
#' @description
#' Calculates empirical probability weighted moments of specific order(s).
#' @param x numeric vector of data.
#' @param order integer, order of probability weighted moment, can be a set of \{0,1,...\}.
#' @param na.rm logical, indicates if NAs should be removed.
#'
#' @return numeric vector, empirical PWM of orders \code{order} of \code{x}.
#' @seealso \code{\link{PWMs}}
#'
#' @examples
#' PWM(rnorm(25))
#' PWM(rnorm(25), order = 2)
#' PWM(rnorm(25), order = c(0, 2, 4))
#' @export
PWM <- function(x, order = 0, na.rm = FALSE) {
  if (is.double(order)) order <- as.integer(order)
  if (!is.integer(order)) stop("order has to be integer type!")
  if (max(order) > length(x)) stop("order higher than number of data points!")
  if (na.rm) {
    x <- na.exclude(x)
  }
  r <- vapply(order, function(r) {
    pwm_C(x, r)
  }, numeric(1))
  setNames(r, paste0("beta", order))
}
