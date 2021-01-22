#' Plot a solution path
#' 
#' @param fit IILasso model
#' @param ... parameters of matlines function
#' 
#' @export
#' 
#' @examples 
#' X <- matrix(c(1,2,3,5,4,7,6,8,9,10), nrow=5, ncol=2)
#' b <- matrix(c(-1,1), nrow=2, ncol=1)
#' e <- matrix(c(0,-0.1,0.1,-0.1,0.1), nrow=5, ncol=1)
#' y <- as.numeric(X %*% b + e)
#' fit <- lasso(X, y)
#' pr <- predict_lasso(fit, X)
#' plot_lasso(fit)
plot_lasso <- function(fit, ...) {
  plot_args <- list(x=log(fit$lambda), y=fit$beta[1,],
                    xlim=range(log(fit$lambda), na.rm=TRUE), 
                    ylim=range(fit$beta, na.rm=TRUE),
                    type="n",
                    xlab="log(lambda)", ylab="coefficients")
  specified_args <- list(...)
  if (length(specified_args)>0) {
    plot_args[names(specified_args)] <- specified_args
  }
  do.call("plot", plot_args)
  # plot(x=log(fit$lambda), y=fit$beta[1,],
  #      xlim=range(log(fit$lambda), na.rm=TRUE), 
  #      ylim=range(fit$beta, na.rm=TRUE),
  #      type="n",
  #      xlab="log(lambda)", ylab="coefficients")
  matlines_args <- list(x=log(fit$lambda), y=t(fit$beta), lty=1)
  if (length(specified_args)>0) {
    matlines_args[names(specified_args)] <- specified_args
  }
  do.call("matlines", matlines_args)
  # matlines(log(fit$lambda), t(fit$beta), lty=1, ...)
}

#' Plot a cross validation error path
#' 
#' @param cv_fit cross validated IILasso model
#' @param ... parameters of 
#' 
#' @export
#' 
#' @examples 
#' X <- matrix(c(1,2,3,5,4,7,6,8,9,10), nrow=5, ncol=2)
#' b <- matrix(c(-1,1), nrow=2, ncol=1)
#' e <- matrix(c(0,-0.1,0.1,-0.1,0.1), nrow=5, ncol=1)
#' y <- as.numeric(X %*% b + e)
#' cv_fit <- cv_lasso(X, y, nfolds=5)
#' fit <- cv_fit$fit
#' pr <- predict_lasso(fit, X, cv_fit$lambda.min)
#' plot_cv_lasso(cv_fit)
plot_cv_lasso <- function(cv_fit, ...) {
  plot_args <- list(log(cv_fit$lambda), cv_fit$cve, type="p", 
                    xlim=range(log(cv_fit$lambda), na.rm=TRUE), 
                    ylim=range(c(cv_fit$cvlo, cv_fit$cvup), na.rm=TRUE),
                    col="red", pch=16,
                    xlab="log(lambda)", ylab="Cross Validation Error")
  specified_args <- list(...)
  if (length(specified_args)>0) {
    plot_args[names(specified_args)] <- specified_args
  }
  do.call("plot", plot_args)
  # plot(log(cv_fit$lambda), cv_fit$cve, type="p", 
  #      xlim=range(log(cv_fit$lambda), na.rm=TRUE), 
  #      ylim=range(c(cv_fit$cvlo, cv_fit$cvup), na.rm=TRUE),
  #      col="red", pch=16,
  #      xlab="log(lambda)", ylab="Cross Validation Error")
  suppressWarnings(arrows(x0=log(cv_fit$lambda), x1=log(cv_fit$lambda),
                          y0=cv_fit$cvlo, y1=cv_fit$cvup,
                          code=3, angle=90, col="gray80", length=.05, ...))
  abline(v=log(cv_fit$lambda[cv_fit$lambda.min.index]),lty=2,lwd=.5)
  abline(v=log(cv_fit$lambda[cv_fit$lambda.1se.index]),lty=3,lwd=.5)
}