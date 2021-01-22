#' @title
#' Summary quantiles
#' @description
#' Calculating and printing of summary statistics to a given quantiles-object.
#'
#' @param object object of quantiles.
#' @param ci.level numeric vector of length 1 giving the confidence level (default is 0.9).
#' @param ... additional arguments submitted to \code{est_cov}.
#'
#' @return A \code{summary.quantiles}-object, a list with dimensions \itemize{
#'  \item \code{q}
#'  \item \code{ci.level}
#'  \item \code{ci}
#'  \item \code{cov}
#' }
#' It is printed with \code{print.summary.quantiles}.
#'
#' @seealso \code{\link{quantiles}}, \code{\link{est_cov}}
#' @examples
#' x <- cbind(rgev(100, shape = .2), rgev(100, shape = .2))
#'
#' q <- quantiles(parameters(TLMoments(x[, 1]), "gev"), c(.9, .95, .99))
#' summary(q)
#' summary(q, select = c(.9, .99))
#'
#' q <- quantiles(parameters(TLMoments(x[, 1], rightrim = 1), "gev"), .95)
#' summary(q)
#'
#' q <- quantiles(parameters(TLMoments(x), "gev"), c(.9, .95, .99))
#' summary(q)
#' summary(q, select = .95)
#'
#' q <- quantiles(as.parameters(loc = 10, scale = 5, shape = .3, distr = "gev"), c(.9, .99))
#' summary(q)
#' summary(q, rightrim = 1, set.n = 250)
#'
#' @method summary quantiles
#' @export
summary.quantiles <- function(object, ci.level = .9, ...) {
  if (length(ci.level) != 1 | !is.numeric(ci.level)) stop("ci must be a numeric vector of length 1. ")
  if (!inherits(object, "quantiles")) stop("First argument has to be of class parameters ")

  UseMethod("summary.quantiles")
}

#' @method summary.quantiles numeric
#' @export
summary.quantiles.numeric <- function(object, ci.level = .9, ...) {
  u <- qnorm(1-(1-ci.level)/2)

  # param ci
  cov <- est_cov(object, ...)
  sel <- paste0("q", names(object)) %in% colnames(cov)
  quan_ci <- object[sel] %-+% (u * sqrt(diag(as.matrix(cov))))

  out <- list(
    q = object,
    ci.level = ci.level,
    ci = cbind(LCL = quan_ci[, 1], quantile = object[sel], UCL = quan_ci[, 2]),
    cov = cov
  )

  class(out) <- c("summary.quantiles")
  out
}

#' @method summary.quantiles matrix
#' @export
summary.quantiles.matrix <- function(object, ci.level = .9, ...) {
  u <- qnorm(1-(1-ci.level)/2)

  # lambda ci
  cov <- est_cov(object, ...)
  sel <- paste0("q", paste(rownames(object), col(object), sep = "_")) %in% colnames(cov)
  quan_ci <- as.vector(object[sel]) %-+% (u * sqrt(diag(as.matrix(cov))))

  out <- list(
    q = object,
    ci.level = ci.level,
    ci = cbind(LCL = quan_ci[, 1], quantile = as.vector(object[sel]), UCL = quan_ci[, 2]),
    cov = cov
  )

  class(out) <- c("summary.quantiles")
  out
}

#' @export
print.summary.quantiles <- function(x, ...) {
  if (any(attr(x$q, "source")$func %in% c("as.PWMs", "as.TLMoments", "as.parameters"))) {
    # Theoretical data
    cat("Confidence intervals based on assumptions: \n")
    cat("\tTrue distribution: ", toupper(attr(x$cov, "distribution")), "\n", sep = "")
    cat("\tTL(", paste0(attr(x$cov, "trimmings"), collapse = ","), ") estimation\n", sep = "")
    cat("\tn = ", attr(x$cov, "n"), "\n", sep = "")
  } else {
    # Empirical data
    ns <- attr(x$q, "source")$n
    cat(length(ns), " data row(s) with n = ", paste0(ns, collapse = ", "), ".\n", sep = "")
    cat("TL(", paste0(attr(x$q, "source")$trim, collapse = ","), ") used to generate ", toupper(attr(x$q, "distr")), " parameters to calculate ", paste0(attr(x$q, "p"), collapse = ", "), " quantile estimates. \n", sep = "")
  }
  cat("\n")
  cat("Approximate ", x$ci.level*100, "% confidence interval of quantiles: \n", sep = "")
  print(x$ci)
  #cat("\n")
  #cat("Covariance matrix of quantile estimates: \n", sep = "")
  #print(x$cov)
}
