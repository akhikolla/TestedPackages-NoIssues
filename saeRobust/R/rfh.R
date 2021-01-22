#' Robust Fay Herriot Model
#'
#' User interface to fit robust Fay-Herriot type models. These models are here
#' framed as linear mixed models. The parameter estimation is robust against
#' outliers. Available models are the standard FH model, a spatial extension, a
#' temporal extension and a spatio-temporal extension.
#'
#' @param formula (formula) a formula specifying the fixed effects part of the
#'   model.
#' @param data (data.frame) a data set.
#' @param samplingVar (character) the name of the variable in \code{data}
#'   containing the sampling variances.
#' @param correlation an optional correlation structure, e.g. \link{corSAR1},
#'   for the random effects part of the model. Default is no correlation, i.e. a
#'   random intercept.
#' @param c (numeric) scalar; a multiplyer constant used in the bias correction.
#'   Default is to make no correction for realisations of direct estimator
#'   within \code{c = 1} times the standard deviation of direct estimator.
#' @param ... arguments passed \link{fitGenericModel}
#'
#' @return
#' A list with the following elements:
#'
#' \itemize{
#' \item \code{call} (language) the call generating the value
#' \item \code{formula} (formula) the formula passed as argument
#' \item \code{samplingVar} (numeric) the vector of sampling variances
#' \item \code{coefficients} (numeric) the vector of regression coefficients
#' \item \code{variance} (numeric) the vector of fitted variance parameter(s)
#' \item \code{iterations} (list) reporting each step in the optimisation
#' \item \code{tol} (numeric) the tolerance level used
#' \item \code{maxIter} (numeric) maximum overall allowed iterations
#' \item \code{maxIterParam} (numeric) maximum allowed iterations for model
#' parameters in each overall iteration
#' \item \code{maxIterRe} (numeric) maximum allowed iterations for fitting the
#' random effects
#' \item \code{k} (numeric) tuning constant in influence function
#' \item \code{K} (numeric) additional tuning constant; often derived from
#' \code{k} to scale down the residual variance
#' \item \code{y} (numeric) the response vector
#' \item \code{x} (Matrix) the design matrix
#' \item \code{re} (numeric) the fitted random effects. Can be c(re1, re2)
#' \item \code{reblup} (numeric) the robust best linear unbiased prediction
#' under the fitted model
#' \item \code{residuals} (numeric) the realised sampling errors
#' \item \code{fitted} (numeric) the fitted values using only the fixed effects
#' part
#' }
#'
#' @details
#'
#' To trigger the spatial and temporal extensions you can supply an argument
#' \code{correlation}. When \code{corSAR1} is used the model of Petrucci and
#' Salvati (2006); for \code{corAR1} the model of Rao and Yu (1994) is used; and
#' for \code{corSAR1AR1} the model of Marhuenda et al. (2013).
#'
#' The methods introducing the robust framework underpinning this implementation
#' can be found in Warnholz (2016). They are based on the results by Sinha and
#' Rao (2009) and Richardson and Welsh (1995).
#'
#' @references
#'
#' Marhuenda, Y., I. Molina and D. Morales (2013). "Small area estimation with
#' spatio-temporal Fay-Herriot models". In: Computational Statistics and Data
#' Analysis 58, pp. 308–325.
#'
#' Pratesi, M. and N. Salvati (2008). "Small area estimation: the EBLUP
#' estimator based on spatially correlated random area effects". In: Statistical
#' Methods & Applications 17, pp. 113–141.
#'
#' Rao, J. N. K. and M. Yu (1994). "Small-Area Estimation by Combining
#' Time-Series and Cross-Sectional Data". In: Canadian Journal of Statistics
#' 22.4, pp. 511–528.
#'
#' Richardson, A. M. and A. H. Welsh (1995). "Robust Restricted Maximum
#' Likelihood in Mixed Linear Models". In: Biometrics 51 (4), pp. 1429–1439.
#'
#' Sinha, S. K. and J. N. K. Rao (2009). "Robust Small Area Estimation". In: The
#' Canadian Journal of Statistics 37 (3), pp. 381–399.
#'
#' Warnholz, S. (2016): "Small Area Estimaiton Using Robust Extension to Area
#' Level Models". Not published (yet).
#'
#' @rdname rfh
#'
#' @export
#' @examples
#'
#' # Non-temporal models:
#' data("grapes", package = "sae")
#' data("grapesprox", package = "sae")
#'
#' fitRFH <- rfh(
#'   grapehect ~ area + workdays - 1,
#'   data = grapes,
#'   samplingVar = "var"
#' )
#'
#' fitRFH
#' summary(fitRFH)
#'
#' plot(fitRFH)
#' plot(predict(fitRFH))
#' plot(mse(fitRFH))
#'
#' \dontrun{
#' # And the same including a spatial structure:
#' fitRSFH <- rfh(
#'   grapehect ~ area + workdays - 1,
#'   data = grapes,
#'   samplingVar = "var",
#'   corSAR1(as.matrix(grapesprox))
#' )
#'
#' # Use the same methods, e.g. plot, for all these implementations:
#' data("spacetime", package = "sae")
#' data("spacetimeprox", package = "sae")
#' nTime <- length(unique(spacetime$Time))
#'
#' fitRTFH <- rfh(
#'   Y ~ X1 + X2,
#'   spacetime,
#'   "Var",
#'   corAR1(nTime = nTime)
#' )
#'
#' fitRSTFH <- rfh(
#'   Y ~ X1 + X2,
#'   spacetime,
#'   "Var",
#'   corSAR1AR1(W = as.matrix(spacetimeprox), nTime = nTime)
#' )
#' }
rfh(formula, data, samplingVar, correlation = NULL, ...) %g% standardGeneric("rfh")

#' @name rfh
#' @usage \S4method{rfh}{formula,data.frame,character,ANY}(formula, data,
#'   samplingVar, correlation, ...)
#' @aliases rfh,formula,data.frame,character,ANY-method
#' @rdname rfh
NULL

rfh(formula ~ formula, data ~ data.frame, samplingVar ~ character, correlation ~ ANY, ...) %m% {
  call <- match.call()
  xy <- makeXY(formula, data)
  samplingVar <- check$samplingVar(data[[samplingVar]])

  stripSelf(retList(
    "rfh",
    public = c("call", "formula"),
    super = rfh(xy$y, xy$x, samplingVar, correlation, ...)
  ))

}

#' @rdname fit
#' @export
rfh(formula ~ numeric, data ~ matrix | Matrix, samplingVar ~ numeric, correlation ~ 'NULL', ...) %m% {
  fitrfh(formula, data, samplingVar, ...)
}

#' @rdname fit
#' @export
rfh(formula ~ numeric, data ~ matrix | Matrix, samplingVar ~ numeric, correlation ~ corSAR1, ...) %m% {
  fitrsfh(formula, data, samplingVar, correlation@W, ...)
}

#' @rdname fit
#' @export
rfh(formula ~ numeric, data ~ matrix | Matrix, samplingVar ~ numeric, correlation ~ corAR1, ...) %m% {
  fitrtfh(formula, data, samplingVar, correlation@nTime, ...)
}

#' @rdname fit
#' @export
rfh(formula ~ numeric, data ~ matrix | Matrix, samplingVar ~ numeric, correlation ~ corSAR1AR1, ...) %m% {
  fitrstfh(formula, data, samplingVar, correlation@W, correlation@nTime, ...)
}

#' @param object (rfh) an object of class rfh
#' @param type (character) one or more in \code{c("linear", "reblup", "reblupbc")}
#'
#' @rdname rfh
#' @export
predict.fitrfh <- function(object, type = "reblup", c = 1, ...) {

  re <- as.numeric(variance(object)$Z() %*% object$re)
  Xb <- fitted.values(object)
  Wbc <- weights(object, c = c)$Wbc

  out <- data.frame(re = re, direct = object$y)
  if (is.element("linear", type)) out$linear <- Xb
  if (is.element("reblup", type)) out$reblup <- Xb + re
  if (is.element("reblupbc", type)) out$reblupbc <- as.numeric(Wbc %*% object$y)

  addAttr(out, c("prediction.fitrfh", "data.frame"), "class")

}
