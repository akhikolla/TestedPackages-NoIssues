#' @title
#' Convert to TLMoments-object
#' @description
#' Convert vector, matrix, list, or data.frame of TL-moments
#' or TL-moment ratios or a PWMs-object to a TLMoments-object in order
#' to be used with TLMoments-functions.
#' The first position of a vector or the first row of a matrix is
#' always used as the L1-moment. The \code{ratios} argument determines if the
#' following positions or rows are used as TL-moments oder TL-moments ratios.
#' The trimming has to be given using the \code{leftrim} and \code{rightrim} arguments.
#'
#' @param x vector or matrix of TL-moments (or TL-moment ratios if ratios is TRUE) or a PWMs-object.
#' The first position or row is always used as the L1-moment, if vector or matrix are given.
#' @param formula if \code{x} is data.frame. See examples.
#' @param ratios boolean, if TRUE the non-first positions or rows of x give
#' L-moment ratios, if FALSE (default) they give L-moments.
#' If ratios are used and the first position or row is NA, L1 is assumed to be 1!
#' @param leftrim,rightrim integer, order of trimmed L-moments.
#' @param ... additional arguments.
#'
#' @return object of class TLMoments, see PWMs help page.
#' @seealso \code{\link{TLMoments}}
#'
#' @examples
#' ### Vector or matrix as input
#' xmat <- cbind(c(1, .2, .05), c(1, .2, NA), c(1.3, NA, .1))
#' xvec <- xmat[, 1]
#' xlist <- lapply(1:ncol(xmat), function(i) xmat[, i])
#' xdat <- data.frame(
#'  station = rep(letters[1:3], each = 1),
#'  season = c("S", "W", "S"),
#'  L1 = c(1, 1, 1.3),
#'  L2 = c(.2, .2, .3),
#'  L3 = c(.05, .04, .1)
#' )
#'
#' as.TLMoments(xvec, rightrim = 1)
#' as.TLMoments(xmat, rightrim = 1)
#' as.TLMoments(xlist, rightrim = 1)
#' as.TLMoments(xdat, cbind(L1, L2, L3) ~ station)
#' as.TLMoments(xdat, .~station+season)
#' as.TLMoments(xdat, cbind(L1, L2, L3) ~ .)
#'
#' parameters(as.TLMoments(xvec, rightrim = 0), "gev")
#' #lmomco::lmom2par(lmomco::vec2lmom(c(1, .2, .25)), "gev")$para
#'
#' xmat <- cbind(c(NA, .2, -.05), c(NA, .2, .2))
#' xvec <- xmat[, 1]
#'
#' as.TLMoments(xvec, ratios = TRUE)
#' as.TLMoments(xmat, ratios = TRUE)
#' parameters(as.TLMoments(xvec, ratios = TRUE), "gev")
#' #lmomco::lmom2par(lmomco::vec2lmom(c(1, .2, -.05)), "gev")$para
#'
#' xmat <- cbind(c(10, .2, -.05), c(10, .2, .2))
#' xvec <- xmat[, 1]
#'
#' as.TLMoments(xvec, ratios = TRUE)
#' as.TLMoments(xmat, ratios = TRUE)
#' parameters(as.TLMoments(xvec, ratios = TRUE), "gev")
#' #lmomco::lmom2par(lmomco::vec2lmom(c(10, .2, -.05)), "gev")$para
#'
#' @rdname as.TLMoments
#' @export
as.TLMoments <- function(x, ..., leftrim, rightrim, ratios) UseMethod("as.TLMoments")

#' @describeIn as.TLMoments as.TLMoments for numeric data vectors
#' @method as.TLMoments numeric
#' @export
as.TLMoments.numeric <- function(x, leftrim = 0L, rightrim = 0L, ratios = FALSE, ...) {
  dim(x) <- c(length(x), 1)
  out <- as.TLMoments(x, ratios = ratios, leftrim = leftrim, rightrim = rightrim)
  out$lambdas <- drop(out$lambdas)
  out$ratios <- drop(out$ratios)
  out
}


#' @describeIn as.TLMoments as.TLMoments for numeric data matrices
#' @method as.TLMoments matrix
#' @export
as.TLMoments.matrix <- function(x, leftrim = 0L, rightrim = 0L, ratios = FALSE, ...) {

  if (!ratios) {
    lambdas <- x
    ratios <- apply(x, 2, calcRatios)
  } else {
    ratios <- x
    lambdas <- apply(x, 2, function(y) calcLambdas(y[-1], L1 = y[1]))
  }

  rownames(lambdas) <- paste0("L", 1L:nrow(lambdas))
  rownames(ratios) <- paste0("T", 1L:nrow(ratios))
  out <- list(
    lambdas = lambdas,
    ratios = ratios
  )

  returnTLMoments(out, leftrim, rightrim, 1L:nrow(lambdas),
                  func = "as.TLMoments",
                  input = x)
}

#' @describeIn as.TLMoments as.TLMoments for numeric data lists
#' @method as.TLMoments list
#' @export
as.TLMoments.list <- function(x, leftrim = 0L, rightrim = 0L, ratios = FALSE, ...) {

  if (!ratios) {
    lambdas <- x
    ratios <- lapply(x, calcRatios)
  } else {
    ratios <- x
    lambdas <- lapply(x, function(y) calcLambdas(y[-1], L1 = y[1]))
  }

  out <- list(
    lambdas = lapply(lambdas, function(x) setNames(x, paste0("L", seq_along(x)))),
    ratios = lapply(ratios, function(x) setNames(x, paste0("T", seq_along(x))))
  )

  returnTLMoments(out, leftrim, rightrim, 1L:length(lambdas[[1]]),
                  func = "as.TLMoments",
                  input = x)
}

#' @describeIn as.TLMoments as.TLMoments for numeric data.frames
#' @method as.TLMoments data.frame
#' @export
as.TLMoments.data.frame <- function(x, formula, leftrim = 0L, rightrim = 0L, ratios = FALSE, ...) {

  nam <- getFormulaSides(formula, names(x))

  if (!ratios) {
    lambdas <- t(as.matrix(x[names(x) %in% nam$lhs]))
    ratios <- apply(lambdas, 2, calcRatios)
  } else {
    ratios <- t(as.matrix(x[names(x) %in% nam$lhs]))
    lambdas <- apply(ratios, 2, function(y) calcLambdas(y[-1], L1 = y[1]))
  }

  out <- list(
    as.data.frame(cbind(x[!(names(x) %in% nam$lhs)], t(lambdas))),
    as.data.frame(cbind(x[!(names(x) %in% nam$lhs)], t(ratios)[, -1]))
  )

  returnTLMoments(out, leftrim, rightrim, 1L:nrow(lambdas),
                  func = "as.TLMoments",
                  input = x,
                  formula = formula)
}
