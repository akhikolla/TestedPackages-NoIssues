#' Fit model on Bootstrap sample
#'
#' These functions help to repeatedly fit a \link{rfh} model on bootstrap
#' samples. Use \code{bootstrap} as a user interface. \code{boot} can be used to
#' extend the framework but is not meant to be used interactively. If you are
#' interested in the parameteric bootstrap for a 'rfh' model you can use the
#' implementation in \link{mse}.
#'
#' @param object a fitted object
#' @param matV the variance of a fitted object used to draw samples. In most
#'   cases this is \code{object}. Alternatively it may be useful to use a
#'   non-robust model.
#' @param B the number of repetitions
#' @param filter a vector indicating which elements in the fittedd object to
#'   keep in each repetition.
#' @param postProcessing a function to process the results. Is applied before
#'   the filter.
#' @param ... arguments passed down to methods
#'
#' @export
#' @rdname bootstrap
#'
#' @examples
#' data(milk, package = "sae")
#' milk$samplingVar <- milk$SD^2
#' modelFit <- rfh(yi ~ as.factor(MajorArea), milk, "samplingVar")
#' bootstrapCoefs <- bootstrap(modelFit, B = 2, filter = "coefficients")
#' do.call(rbind, unlist(bootstrapCoefs, FALSE))
bootstrap <- function(object, matV = variance(object), B = NULL, ...) {
  # This function avoids the 'unecessary' methods for missing values when we
  # want S4 dispatch with default values in the generic
  boot(object, matV, B, ...)
}

#' @export
#' @rdname bootstrap
boot(object, matV, B, ...) %g% standardGeneric("boot")

#' @export
#' @rdname bootstrap
boot(object, matV, B ~ integer|numeric, filter = NULL, postProcessing = identity, ...) %m% {
  if (is.null(filter)) {
    pbreplicate(B, postProcessing(boot(object, matV, NULL, ...)), FALSE)
  } else {
    pbreplicate(B, postProcessing(boot(object, matV, NULL, ...))[filter], FALSE)
  }
}

#' @export
#' @rdname bootstrap
boot(object ~ rfh, matV ~ rfhVariance, B ~ NULL, ...) %m% {
  # This we need to get directly in the update method for fitrfh class
  # Otherwise weired things are happening in the call for the S4-generic
  class(object) <- class(object)[-1] # this hopefully only removes 'rfh'

  # Bootstrap sample:
  Xb <- fitted.values(object)
  re <- MASS::mvrnorm(1, mu = rep(0, nrow(matV$Vu())), matV$Vu())
  e <- MASS::mvrnorm(1, mu = rep(0, length(Xb)), matV$Ve())
  trueY <- as.numeric(Xb + matV$Z() %*% re)
  y <- trueY + e

  # refit:
  out <- update(object, y = y)
  out$trueY <- trueY
  out
}

