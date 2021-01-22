#' @title
#' Calculation of regionalized TL-moments
#' @description
#' regionalize takes the result of TLMoments and calculates a weighted mean
#' of TL-moments and TL-moment ratios.
#'
#' @param x object returned by TLMoments.
#' @param w numeric vector giving the weights. Default: Sample lengths of corresponding
#' data. Internally scaled so that it adds up to 1.
#' @param reg.lambdas logical, if TRUE (default) regionalization is
#' based upon TL-moments. If false it's based on TL-moment-ratios.
#' @param ... additional arguments, not used at the moment.
#'
#' @return list of two dimensions: \code{lambdas}/\code{ratios} are numeric vectors
#' consisting of the regionalized TL-moments/TL-moment-ratios. The list has
#' the class \code{TLMoments}. The object contains the following attributes: \itemize{
#'  \item \code{leftrim}: a numeric giving the used leftrim-argument
#'  \item \code{rightrim}: a numeric giving the used rightrim-argument
#'  \item \code{order}: a integer vector with corresponding TL-moment orders
#'  \item \code{source}: a list with background information (used function, data,
#'  n, formula, computation.method; mainly for internal purposes)
#' }
#'
#' @examples
#' xmat <- matrix(rgev(100), nc = 4)
#' xvec <- xmat[, 3]
#' xlist <- lapply(1L:ncol(xmat), function(i) xmat[, i])
#' xdat <- data.frame(
#'  station = rep(letters[1:2], each = 50),
#'  season = rep(c("S", "W"), 50),
#'  hq = as.vector(xmat)
#' )
#'
#' regionalize(TLMoments(xmat))
#' regionalize(TLMoments(xlist))
#' regionalize(TLMoments(xdat, hq ~ station))
#' # For numeric vector TLMoments, nothing happens:
#' regionalize(TLMoments(xvec))
#'
#' tlm <- TLMoments(xmat)
#' regionalize(tlm)
#' regionalize(tlm, reg.lambdas = FALSE)
#'
#' parameters(regionalize(tlm), "gev")
#' parameters(regionalize(tlm, reg.lambdas = FALSE), "gev")
#'
#' quantiles(parameters(regionalize(tlm), "gev"), c(.99, .999))
#' quantiles(parameters(regionalize(tlm, reg.lambdas = FALSE), "gev"), c(.99, .999))
#'
#'
#' # With magrittr
#' library(magrittr)
#' matrix(rgev(200, shape = .3), nc = 5) %>%
#'  TLMoments(rightrim = 1) %>%
#'  regionalize %>%
#'  parameters("gev") %>%
#'  quantiles(c(.99, .999))
#' @rdname regionalize
#' @export
regionalize <- function(x, ...)  {
  if (!inherits(x, "TLMoments"))
    stop("x must be object of class TLMoments. ")

  UseMethod("regionalize", x$lambdas)
}

#' @rdname regionalize
#' @method regionalize numeric
#' @export
regionalize.numeric <- function(x, ...) {
  x
}

#' @rdname regionalize
#' @method regionalize matrix
#' @export
regionalize.matrix <- function(x, w = attr(x, "source")$n, reg.lambdas = TRUE, ...) {

  # Ensure that the sum of weights is 1.
  if (sum(w) != 1) w <- w / sum(w)

  # Two versions:
  # 1) Regionalize lambdas and calculate taus (default)
  # 2) Regionalize taus (and l1) und calculate lambdas
  if (reg.lambdas) {
    lambdas <- apply(x$lambdas, 1, stats::weighted.mean, w = w)
    ratios <- calcRatios(lambdas)
  } else {
    ratios <- apply(x$ratios, 1, stats::weighted.mean, w = w)
    l1 <- weighted.mean(x$lambdas[1, ], w)
    lambdas <- calcLambdas(ratios[-1], l1)
 }
  out <- list(lambdas = lambdas, ratios = ratios)

  do.call(returnTLMoments, c(
    list(out = out, leftrim = attr(x, "leftrim"),
         rightrim = attr(x, "rightrim"), order = seq_along(lambdas)),
    func = "regionalize",
    reg.weights = list(w),
    attr(x, "source")
  ))
}

#' @rdname regionalize
#' @method regionalize data.frame
#' @export
regionalize.data.frame <- function(x, w = attr(x, "source")$n, reg.lambdas = TRUE, ...) {

  # Ensure that the sum of weights is 1.
  if (sum(w) != 1) w <- w / sum(w)

  # Two versions:
  # 1) Regionalize lambdas and calculate taus (default)
  # 2) Regionalize taus (and l1) und calculate lambdas
  nam <- getFormulaSides(attr(x, "source")$formula)
  if (reg.lambdas) {
    lambdas <- t(x$lambdas[!(names(x$lambdas) %in% nam$rhs)])
    lambdas <- apply(lambdas, 1, stats::weighted.mean, w = w)
    ratios <- calcRatios(lambdas)
  } else {
    ratios <- t(x$ratios[!(names(x$ratios) %in% nam$rhs)])
    ratios <- apply(ratios, 1, stats::weighted.mean, w = w)
    l1 <- weighted.mean(x$lambdas[!(names(x$lambdas) %in% nam$rhs)][, 1], w)
    lambdas <- calcLambdas(ratios, l1)
  }
  out <- list(lambdas = lambdas, ratios = ratios)

  do.call(returnTLMoments, c(
    list(out = out, leftrim = attr(x, "leftrim"),
         rightrim = attr(x, "rightrim"), order = seq_along(lambdas)),
    func = "regionalize",
    reg.weights = list(w),
    attr(x, "source")
  ))
}

#' @rdname regionalize
#' @method regionalize list
#' @export
regionalize.list <- function(x, w = attr(x, "source")$n, reg.lambdas = TRUE, ...) {

  # Ensure that the sum of weights is 1.
  if (sum(w) != 1) w <- w / sum(w)

  # Two versions:
  # 1) Regionalize lambdas and calculate taus (default)
  # 2) Regionalize taus (and l1) und calculate lambdas
  if (reg.lambdas) {
    lambdas <- do.call(cbind, x$lambdas)
    lambdas <- apply(lambdas, 1, stats::weighted.mean, w = w)
    ratios <- calcRatios(lambdas)
  } else {
    ratios <- do.call(cbind, x$ratios)
    ratios <- apply(ratios, 1, stats::weighted.mean, w = w)
    l1 <- weighted.mean(sapply(x$lambdas, getElement, 1), w)
    lambdas <- calcLambdas(ratios[-1], l1)
  }
  out <- list(lambdas = lambdas, ratios = ratios)

  do.call(returnTLMoments, c(
    list(out = out, leftrim = attr(x, "leftrim"),
         rightrim = attr(x, "rightrim"), order = seq_along(lambdas)),
    func = "regionalize",
    reg.weights = list(w),
    attr(x, "source")
  ))
}
