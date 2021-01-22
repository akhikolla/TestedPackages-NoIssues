#' @title
#' Summary parameters
#' @description
#' Calculating and printing of summary statistics to a given parameters-object.
#'
#' @param object object of parameters.
#' @param ci.level numeric vector of length 1 giving the confidence level (default is 0.9).
#' @param ... additional arguments submitted to \code{est_cov}.
#' @return A \code{summary.parameters}-object, a list with dimensions \itemize{
#'  \item \code{param}
#'  \item \code{ci.level}
#'  \item \code{ci}
#'  \item \code{cov}
#' }
#' It is printed with \code{print.summary.parameters}.
#'
#' @seealso \code{\link{parameters}}, \code{\link{est_cov}}
#' @examples
#' x <- cbind(rgev(100, shape = .2), rgev(100, shape = .2))
#'
#' p <- parameters(TLMoments(x[, 1]), "gev")
#' summary(p)
#' summary(p, select = c("scale", "shape"))
#'
#' p <- parameters(TLMoments(x[, 1], rightrim = 1), "gev")
#' summary(p)
#'
#' p <- parameters(TLMoments(x), "gev")
#' summary(p)
#' summary(p, select = "shape")
#'
#' p <- as.parameters(loc = 10, scale = 5, shape = .3, distr = "gev")
#' summary(p)
#' summary(p, rightrim = 1, set.n = 250)
#'
#' @method summary parameters
#' @export
summary.parameters <- function(object, ci.level = .9, ...) {
  if (length(ci.level) != 1 | !is.numeric(ci.level)) stop("ci must be a numeric vector of length 1. ")
  if (!inherits(object, "parameters")) stop("First argument has to be an object of class parameters. ")

  UseMethod("summary.parameters")
}

#' @method summary.parameters numeric
#' @export
summary.parameters.numeric <- function(object, ci.level = .9, ...) {
  u <- qnorm(1-(1-ci.level)/2)

  # param ci
  cov <- est_cov(object, ...)
  sel <- names(object) %in% colnames(cov)
  param_ci <- object[sel] %-+% (u * sqrt(diag(as.matrix(cov))))

  out <- list(
    param = object,
    ci.level = ci.level,
    ci = cbind(LCL = param_ci[, 1], param = object[sel], UCL = param_ci[, 2]),
    cov = cov
  )

  class(out) <- c("summary.parameters")
  out
}

#' @method summary.parameters matrix
#' @export
summary.parameters.matrix <- function(object, ci.level = .9, ...) {
  u <- qnorm(1-(1-ci.level)/2)

  # lambda ci
  cov <- est_cov(object, ...)
  sel <- paste(rownames(object), col(object), sep = "_") %in% colnames(cov)
  param_ci <- as.vector(object[sel]) %-+% (u * sqrt(diag(as.matrix(cov))))

  out <- list(
    param = object,
    ci.level = ci.level,
    ci = cbind(LCL = param_ci[, 1], param = as.vector(object[sel]), UCL = param_ci[, 2]),
    cov = cov
  )

  class(out) <- c("summary.parameters")
  out
}

#' @export
print.summary.parameters <- function(x, ...) {
  # summary prints an overview of the data and what was done and ci of the parameters

  if (any(attr(x$param, "source")$func %in% c("as.PWMs", "as.TLMoments", "as.parameters"))) {
    # Theoretical data
    cat("Distribution: ", toupper(attr(x$cov, "distribution")), "(", paste(names(x$param), round(x$param, 4), sep = "=", collapse = ", "), ")\n", sep = "")
    cat("Confidence intervals based on assumptions: \n")
    cat("\tTL(", paste0(attr(x$cov, "trimmings"), collapse = ","), ")\n", sep = "")
    cat("\tn = ", attr(x$cov, "n"), "\n", sep = "")
  } else {
    # Empirical data
    ns <- attr(x$param, "source")$n
    cat(length(ns), " data row(s) with n = ", paste0(ns, collapse = ", "), ".\n", sep = "")
    cat("TL(", paste0(attr(x$param, "source")$trim, collapse = ","), ") used to generate ", toupper(attr(x$param, "distribution")), " parameters. \n", sep = "")
  }
  cat("\n")
  cat("Approximate ", x$ci.level*100, "% confidence interval of parameters: \n", sep = "")
  print(x$ci)
}
