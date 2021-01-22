#' Compute the Mean Squared Error of an Estimator
#'
#' A generic function to compute the mean squared error of the predicted values
#' under the estimated model. See also \link{rfh} for examples.
#'
#' @param object (see methods) an object containing the estimation result, e.g.
#'   \link{rfh}
#' @param type (character) the type of the MSE. Available are 'pseudo' and
#'   'boot'
#' @param B (numeric) number of bootstrap repetitions
#' @param predType (character) the type of prediction: \code{c("reblup",
#'   "reblupbc")}
#' @param ... arguments passed to methods
#'
#' @details
#'
#' Type pseudo is an approximation of the MSE based on a pseudo
#' linearisation approach by Chambers, et. al. (2011). The specifics can be
#' found in Warnholz (2016). Type boot implements a parameteric bootstrap for
#' these methods.
#'
#' @references
#'
#' Chambers, R., H. Chandra and N. Tzavidis (2011). "On bias-robust mean squared
#' error estimation for pseudo-linear small area estimators". In: Survey
#' Methodology 37 (2), pp. 153–170.
#'
#' Warnholz, S. (2016): "Small Area Estimaiton Using Robust Extension to Area
#' Level Models". Not published (yet).
#'
#' @export
#' @rdname mse
#'
#' @examples
#' data("grapes", package = "sae")
#' data("grapesprox", package = "sae")
#'
#' fitRFH <- rfh(
#'   grapehect ~ area + workdays - 1,
#'   data = grapes,
#'   samplingVar = "var"
#' )
#'
#' mseRFH <- mse(fitRFH)
#' plot(mseRFH)
mse <- function(object, ...) UseMethod("mse")

#' @export
#' @rdname mse
mse.fitrfh <- function(object, type = "pseudo", predType = "reblupbc", B = 100, ...) {

  pseudo <- function(Xb, matV, Wtype = "W") {
    # Compute:
    # Var: A Z Vu t(Z) t(A) + W Ve t(W)
    # Bias: W Xb - Xb
    # diag(var) + bias^2
    W <- weights(object)[[Wtype]]
    G <- matV$Z() %*% tcrossprod(matV$Vu(), matV$Z())
    A <- W - Diagonal(length(object$y))
    var <- A %*% tcrossprod(G, A) + W %*% tcrossprod(matV$Ve(), W)
    bias <- W %*% Xb - Xb
    as.numeric(diag(var) + bias^2)
  }

  boot <- function(B, matV) {

    resList <- bootstrap(
      object, matV, B = B,
      postProcessing = function(fit) c(predict(fit, type = c("reblup", "reblupbc")), fit["trueY"])
    )

    reblup <- colMeans(do.call(rbind, lapply(resList, function(el) {
      (el$reblup - el$trueY)^2
    })))

    reblupbc <- colMeans(do.call(rbind, lapply(resList, function(el) {
      (el$reblupbc - el$trueY)^2
    })))

    retList(public = c("reblup", "reblupbc"))

  }

  matV <- variance(object)
  Xb <- fitted.values(object)

  out <- predict(object, type = predType)[c("direct", predType)]
  out$samplingVar <- object$samplingVar
  if ("pseudo" %in% type && "reblup" %in% predType) out$pseudo <- pseudo(Xb, matV)
  if ("pseudo" %in% type && "reblupbc" %in% predType) out$pseudobc <- pseudo(Xb, matV, "Wbc")
  if ("boot" %in% type) {
    boots <- boot(B, matV)
    if ("reblup" %in% predType) out$boot <- boots$reblup
    if ("reblupbc" %in% predType) out$bootbc <- boots$reblupbc
  }

  addAttr(out, c("mse.fitrfh", "data.frame"), "class")

}
