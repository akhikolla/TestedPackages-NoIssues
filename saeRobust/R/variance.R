#' Construct variance
#'
#' A generic function to construct the different variance components of an
#' object. You may want to use this in conjunction with \link{bootstrap}.
#'
#' @param .object,object an object
#' @param c (numeric) scalar
#' @param ... arguments passed to method
#'
#' @export
#' @rdname variance
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
#' # The variance component of a mixed linear model:
#' matV <- variance(fitRFH)
#' # The full variance matrix:
#' matV$V()
#'
#' # The sampling error component
#' matV$Ve()
#'
#' # the random effects component
#' matV$Vu()
variance <- function(.object, ...) UseMethod("variance")

#' @export
#' @rdname variance
variance.fitrfh <- function(.object, ...) {

  expose(matVFH(.object$variance, .object$samplingVar))

  retList("rfhVariance")
}

#' @export
#' @rdname variance
variance.fitrsfh <- function(.object, ...) {

  expose(matVSFH(
    .object$variance[1],
    .object$variance[2],
    .object$W,
    .object$samplingVar))

  retList("rfhVariance")
}

#' @export
#' @rdname variance
variance.fitrtfh <- function(.object, ...) {

  expose(matVTFH(
    .object$variance[1],
    .object$variance[c(2, 3)],
    .object$nTime,
    .object$samplingVar))

  retList("rfhVariance")
}

#' @export
#' @rdname variance
variance.fitrstfh <- function(.object, ...) {

  expose(matVSTFH(
    .object$variance[1:2],
    .object$variance[3:4],
    .object$W,
    .object$nTime,
    .object$samplingVar))

  retList("rfhVariance")

}

#' @export
#' @rdname variance
weights.fitrfh <- function(object, c = 1, ...) {

  V <- variance(object)

  W <-  matW(
    y = object$y,
    X = object$x,
    beta = object$coefficients,
    re = object$re,
    matV = V,
    psi = . %>% psiOne(object$k)
  )

  A <- matA(
    y = object$y,
    X = object$x,
    beta = object$coefficients,
    matV = V,
    psi = . %>% psiOne(object$k)
  )

  B <- matB(
    y = object$y,
    X = object$x,
    beta = object$coefficients,
    re = object$re,
    V,
    psi = . %>% psiOne(object$k)
  )

  Wbc <- matWbc(
    y = object$y,
    reblup = object$reblup,
    W = W,
    samplingVar = object$samplingVar,
    c = c
  )

  stripSelf(retList("rfhWeights", c("W", "Wbc", "A", "B")))

}

