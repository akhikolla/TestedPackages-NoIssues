#' Plot a solution path
#'
#' @param x hmlasso model
#' @param xlim x range
#' @param ylim y range
#' @param ... parameters of matlines function
#'
#' @examples
#' X_incompl <- as.matrix(iris[, 1:3])
#' X_incompl[1:5,1] <- NA
#' X_incompl[6:10,2] <- NA
#' y <- iris[, 4]
#' fit <- hmlasso(X_incompl, y, nlambda=50, lambda.min.ratio=1e-2)
#' plot(fit)
#'
#' @export
plot.hmlasso <- function(x, xlim=NULL, ylim=NULL, ...) {
  if (all(is.na(x$beta) | is.infinite(x$beta))) {
    warning("beta is all NA or infinite")
    return()
  }
  if (missing(xlim) & missing(ylim)) {
    plot(x=log(x$lambda), y=x$beta[1,],
         xlim=range(log(x$lambda), na.rm=TRUE, finite=TRUE),
         ylim=range(x$beta, na.rm=TRUE, finite=TRUE),
         type="n",
         xlab="log(lambda)", ylab="coefficients")
  } else {
    plot(x=log(x$lambda), y=x$beta[1,],
         xlim=xlim, ylim=ylim,
         type="n",
         xlab="log(lambda)", ylab="coefficients")
  }
  matlines(log(x$lambda), t(x$beta), lty=1, ...)
}

#' Plot a cross validation error path
#'
#' @param x cross validated hmlasso model
#' @param xlim x range
#' @param ylim y range
#' @param ... parameters of
#'
#' @examples
#' X_incompl <- as.matrix(iris[, 1:3])
#' X_incompl[1:5,1] <- NA
#' X_incompl[6:10,2] <- NA
#' y <- iris[, 4]
#' cv_fit <- cv.hmlasso(X_incompl, y, nlambda=50, lambda.min.ratio=1e-2)
#' plot(cv_fit)
#' plot(cv_fit$fit)
#'
#' @export
plot.cv.hmlasso <- function(x, xlim=NULL, ylim=NULL, ...) {
  if (all(is.na(x$cve) | is.infinite(x$cve))) {
    warning("beta is all NA or infinite")
    return()
  }
  if (missing(xlim) & missing(ylim)) {
    plot(log(x$lambda), x$cve, type="p",
         xlim=range(log(x$lambda), na.rm=TRUE, finite=TRUE),
         ylim=range(c(x$cvlo, x$cvup), na.rm=TRUE, finite=TRUE),
         col="red", pch=16,
         xlab="log(lambda)", ylab="Cross Validation Error")
  } else {
    plot(log(x$lambda), x$cve, type="p",
         xlim=xlim, ylim=ylim,
         col="red", pch=16,
         xlab="log(lambda)", ylab="Cross Validation Error")
  }
  suppressWarnings(arrows(x0=log(x$lambda), x1=log(x$lambda),
                          y0=x$cvlo, y1=x$cvup,
                          code=3, angle=90, col="gray80", length=.05, ...))
  abline(v=log(x$lambda[x$lambda.min.index]),lty=2,lwd=.5)
  abline(v=log(x$lambda[x$lambda.1se.index]),lty=3,lwd=.5)
}
