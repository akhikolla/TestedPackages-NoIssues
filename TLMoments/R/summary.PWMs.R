#' @title
#' Summary PWMs
#' @description
#' Calculating and printing of summary statistics to a given PWMs-object.
#'
#' @param object object of PWMs.
#' @param ci.level numeric vector of length 1 giving the confidence level (default is 0.9).
#' @param ... additional arguments submitted to \code{est_cov}.
#'
#' @return A \code{summary.PWMs}-object, a list with dimensions \itemize{
#'  \item \code{pwm}
#'  \item \code{ci.level}
#'  \item \code{ci}
#'  \item \code{cov}
#' }
#' It is printed with \code{print.summary.PWMs}.
#'
#' @seealso \code{\link{PWMs}}, \code{\link{est_cov}}
#' @examples
#' x <- cbind(rgev(100, shape = .2), rgev(100, shape = .2))
#'
#' summary(PWMs(x[, 1]))
#' summary(PWMs(x[, 1]), distr = "gev")
#' summary(PWMs(x[, 1]), distr = "gev", select = 1:2)
#'
#' summary(PWMs(x))
#' summary(PWMs(x), select = 1:2)
#'
#' \dontrun{
#' summary(as.PWMs(c(15, 4, .5)))
#' }
#' @method summary PWMs
#' @export
summary.PWMs <- function(object, ci.level = .9, ...) {
  if (length(ci.level) != 1 | !is.numeric(ci.level))
    stop("ci must be a numeric vector of length 1. ")
  if (!inherits(object, "PWMs"))
    stop("First argument has to be of class parameters ")
  if (any(attr(object, "source")$func %in% c("as.PWMs", "as.TLMoments", "as.parameters")))
    stop("No summary for theoretical PWMs available. ")

  UseMethod("summary.PWMs")
}

#' @method summary.PWMs numeric
#' @export
summary.PWMs.numeric <- function(object, ci.level = .9, ...) {
  u <- qnorm(1-(1-ci.level)/2)

  # param ci
  cov <- est_cov(object, ...)
  sel <- names(object) %in% colnames(cov)
  pwm_ci <- object[sel] %-+% (u * sqrt(diag(as.matrix(cov))))

  out <- list(
    pwm = object,
    ci.level = ci.level,
    ci = cbind(LCL = pwm_ci[, 1], pwm = object[sel], UCL = pwm_ci[, 2]),
    cov = cov
  )

  class(out) <- c("summary.PWMs")
  out
}

#' @method summary.PWMs matrix
#' @export
summary.PWMs.matrix <- function(object, ci.level = .9, ...) {
  u <- qnorm(1-(1-ci.level)/2)

  # lambda ci
  cov <- est_cov(object, ...)
  sel <- paste(rownames(object), col(object), sep = "_") %in% colnames(cov)
  pwm_ci <- as.vector(object[sel]) %-+% (u * sqrt(diag(as.matrix(cov))))

  out <- list(
    pwm = object,
    ci.level = ci.level,
    ci = cbind(LCL = pwm_ci[, 1], pwm = as.vector(object[sel]), UCL = pwm_ci[, 2]),
    cov = cov
  )

  class(out) <- c("summary.PWMs")
  out
}

#' @export
print.summary.PWMs <- function(x, ...) {
  ns <- attr(x$pwm, "source")$n
  cat(length(ns), " data row(s) with n = ", paste0(ns, collapse = ", "), ".\n", sep = "")
  cat("\n")
  cat("Approximate ", x$ci.level*100, "% confidence interval of PWMs: \n", sep = "")
  print(x$ci)
  #cat("\n")
  #cat("Covariance matrix of PWM estimates: \n", sep = "")
  #print(x$cov)
}
