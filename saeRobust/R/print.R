#' @method print rfh
#' @export
print.rfh <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  # Call:

  callChar <- unlist(strsplit(deparse(x$call), split = " "))
  callChar <- callChar[!(callChar == "")]
  cat("\nCall:\n")
  cat(callChar, "\n\n", sep = " ", fill = TRUE)

  # Coefficients:
  if (length(coef(x))) {
    cat("Coefficients:\n")
    printDefault(coefficients(x), digits)
  }
  else cat("No coefficients\n")
  cat("\n")

  # Variance Components:
  cat("Variance Components:\n")
  printDefault(x$variance, digits)
  cat("\n")

  invisible(x)

}

printDefault <- function(x, digits) {
  print.default(
    format(x, digits = digits),
    print.gap = 2L, quote = FALSE
  )
}

#' @method summary rfh
#' @export
summary.rfh <- function(object, ...) {

  # The scores:
  object$score <- score(object)
  names(object$score) <- c("coefficients", "variance", "re")
  names(object$score$coefficients) <- names(coefficients(object))
  names(object$score$variance) <- names(object$variance)

  # Iterations:
  object$iter <- c(
    "model parameter" = NROW(object$iterations$coefficients),
    "random effects" = NROW(object$iterations$re)
  )

  retList("summary.rfh", super = object)
}

#' @method print summary.rfh
#' @export
print.summary.rfh <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  tmp <- x
  class(tmp) <- class(tmp)[-1]
  print(tmp)

  printDefault(digits = digits, rbind(
    "Random Effects" = summary(x$re),
    "Residuals" = summary(residuals(x))
  ))

  cat("\n\n## Solutions to the Robust Estimation Equations:\n\n")
  printDefault(c(x$score$coefficients, x$score$variance), digits)

  cat("\nRandom Effects:\n")
  printDefault(summary(x$score$re), digits)


  cat("\n\n## Iterations:\n")
  cat("Allowed: ", x$maxIter, " (", x$maxIterParam, ") -- ", x$maxIterRe, "\n", sep = "")
  printDefault(x$iter, digits)

  invisible(x)

}
