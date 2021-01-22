#' @title
#' Summary TLMoments
#' @description
#' Calculating and printing of summary statistics to a given TLMoments-object.
#'
#' @param object object of TLMoments.
#' @param ci.level numeric vector of length 1 giving the confidence level (default is 0.9).
#' @param ... additional arguments submitted to \code{est_cov}.
#' @return A \code{summary.TLMoments}-object, a list with dimensions \itemize{
#'  \item \code{tlm}
#'  \item \code{ci.level}
#'  \item \code{lambda.ci}
#'  \item \code{lambda.cov}
#'  \item \code{ratio.ci}
#'  \item \code{ratio.cov}
#' }
#' It is printed with \code{print.summary.TLMoments}.
#'
#' @seealso \code{\link{TLMoments}}, \code{\link{est_cov}}
#' @examples
#' tlm <- TLMoments(rgev(100, shape = .2))
#' summary(tlm)
#'
#' tlm <- TLMoments(rgev(100, shape = .2), rightrim = 1)
#' summary(tlm, select = 3:4)
#'
#' tlm <- TLMoments(rgev(100, shape = .2), max.order = 2, rightrim = 1)
#' summary(tlm)
#'
#' tlm <- TLMoments(matrix(rgev(100, shape = .2), nc = 2))
#' summary(tlm, select = 3:4)
#'
#' tlm <- TLMoments(matrix(rgev(100, shape = .2), nc = 2), max.order = 3)
#' summary(tlm, ci = .95, distr = "gev")
#'
#' tlm <- as.TLMoments(c(15, 5, 1.3))
#' summary(tlm, distr = "gev", set.n = 100)
#'
#' @method summary TLMoments
#' @export
summary.TLMoments <- function(object, ci.level = .9, ...) {
  if (length(ci.level) != 1 | !is.numeric(ci.level)) stop("ci must be a numeric vector of length 1. ")
  if (!inherits(object, "TLMoments")) stop("First argument has to be of class TLMoments. ")

  UseMethod("summary.TLMoments", object$lambdas)
}

#' @method summary.TLMoments numeric
#' @export
summary.TLMoments.numeric <- function(object, ci.level = .9, ...) {
  u <- qnorm(1-(1-ci.level)/2)

  # covs
  cov <- est_cov(object, ...)
  sel <- names(object$lambdas) %in% colnames(cov$lambdas)

  # lambda ci
  lambda_ci <- object$lambdas[sel] %-+% (u * sqrt(diag(cov$lambdas)))
  lambda_ci <- cbind(lambda_ci[, 1], object$lambdas[sel], lambda_ci[, 2])
  colnames(lambda_ci) <- c("LCL", "lambda_hat", "UCL")

  out <- list(
    tlm = object,
    cov = cov,
    ci.level = ci.level,
    lambda.ci = lambda_ci
  )

  # tau ci
  if (length(attr(object, "order")) >= 2) {
    ratio_ci <- object$ratios[sel][-1] %-+% (u * sqrt(diag(as.matrix(cov$ratios))))
    out$ratio.ci <- cbind(ratio_ci[, 1], object$ratios[sel][-1], ratio_ci[, 2])
    colnames(out$ratio.ci) <- c("LCL", "tau_hat", "UCL")
  }

  class(out) <- c("summary.TLMoments")
  out
}

#' @method summary.TLMoments matrix
#' @export
summary.TLMoments.matrix <- function(object, ci.level = .9, ...) {
  u <- qnorm(1-(1-ci.level)/2)

  # covs
  cov <- est_cov(object, ...)
  sel <- paste(rownames(object$lambdas), col(object$lambdas), sep = "_") %in% colnames(cov$lambdas)

  # lambda ci
  lambda_ci <- as.vector(object$lambdas[sel]) %-+% (u * sqrt(diag(cov$lambdas)))
  lambda_ci <- cbind(lambda_ci[, 1], as.vector(object$lambdas[sel]), lambda_ci[, 2])
  colnames(lambda_ci) <- c("LCL", "lambda_hat", "UCL")

  out <- list(
    tlm = object,
    cov = cov,
    ci.level = ci.level,
    lambda.ci = lambda_ci
  )

  # tau ci
  if (length(attr(object, "order")) >= 2) {
    ratio_ci <- as.numeric(na.exclude(object$ratios[sel])) %-+% (u * sqrt(diag(cov$ratios)))
    out$ratio.ci <- cbind(ratio_ci[, 1], as.numeric(na.exclude(object$ratios[sel])), ratio_ci[, 2])
    colnames(out$ratio.ci) <- c("LCL", "tau_hat", "UCL")
  }

  class(out) <- c("summary.TLMoments")
  out
}

#' @export
print.summary.TLMoments <- function(x, ...) {
  if (any(attr(x$tlm, "source")$func %in% c("as.PWMs", "as.TLMoments", "as.parameters"))) {
    # Theoretical data
    cat("TL(", attr(x$tlm, "leftrim"), ",", attr(x$tlm, "rightrim"), ") given.  \n", sep = "")
    cat("Confidence intervals based on assumptions: \n")
    cat("\tTrue distribution: ", toupper(attr(x$cov, "distribution")), "\n", sep = "")
    cat("\tn = ", attr(x$cov, "n"), "\n", sep = "")
  } else {
    # Empirical data
    ns <- attr(x$tlm, "source")$n
    cat(length(ns), " data row(s) with n = ", paste0(ns, collapse = ", "), ".\n", sep = "")
    cat("TL(", attr(x$tlm, "leftrim"), ",", attr(x$tlm, "rightrim"), ") calculated. \n", sep = "")
  }
  cat("\n")
  cat("Approximate ", x$ci.level*100, "% confidence interval of TL moments: \n", sep = "")
  print(x$lambda.ci)
  if (!is.null(x$ratio.ci)) {
    cat("Approximate ", x$ci.level*100, "% confidence interval of TL moment ratios: \n", sep = "")
    print(x$ratio.ci)
  }
}
