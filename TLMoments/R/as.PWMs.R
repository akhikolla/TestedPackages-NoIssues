#' @title
#' Convert to PWMs-object
#' @description
#' Convert vector, matrix, list, or data.frame to PWMs-objects.
#'
#' @param x vector or matrix of PWMs.
#' @param formula if \code{x} is data.frame. See examples.
#' @param order integer, corresponding order to given PWMs. If NULL, order is set to 0:(length(x)-1).
#' @param ... additional arguments.
#'
#' @return object of class PWMs, see PWMs help page.
#' @seealso \code{\link{PWMs}}
#'
#' @examples
#' xmat <- cbind(c(0.12, .41, .38, .33), c(.05, 0.28, .25, .22))
#' xvec <- xmat[, 1]
#' xlist <- lapply(1:ncol(xmat), function(i) xmat[, i])
#' xdat <- data.frame(
#'  station = letters[1:3],
#'  season = c("S", "W", "S"),
#'  b0 = c(.12, .15, .05),
#'  b1 = c(.41, .33, .28),
#'  b2 = c(.38, .18, .25)
#' )
#'
#' as.PWMs(xvec)
#' as.PWMs(xvec[-2], order = c(0, 2, 3))
#'
#' as.PWMs(xmat)
#' as.PWMs(xmat[-2, ], order = c(0, 2, 3))
#'
#' as.PWMs(xlist)
#'
#' as.PWMs(xdat, cbind(b0, b1, b2) ~ station)
#' as.PWMs(xdat, . ~ station + season)
#' as.PWMs(xdat, cbind(b0, b2) ~ station, order = c(0, 2))
#'
#' p <- as.PWMs(xdat, cbind(b0, b1, b2) ~ station)
#' TLMoments(p)
#'
#' (p <- as.PWMs(xdat, cbind(b0, b1) ~ station))
#' #parameters(TLMoments(p), "gev") # => error
#' #parameters(TLMoments(p), "gpd") # => error
#' parameters(TLMoments(p), "gpd", u = 10)
#'
#' (p <- as.PWMs(xdat, cbind(b0, b2) ~ station, order = c(0, 2)))
#' #TLMoments(p) # => error
#'
#' @rdname as.PWMs
#' @export
as.PWMs <- function(x, ..., order = NULL) {
  # TODO: General error catches
  if (!is.null(order)) {
    if (is.double(order)) order <- as.integer(order)
    if (!is.integer(order)) stop("order has to be integer type!")
  }

  UseMethod("as.PWMs")
}


#' @describeIn as.PWMs as.PWMs for numeric data vectors
#' @method as.PWMs numeric
#' @export
as.PWMs.numeric <- function(x, order = seq_along(x)-1, ...) {
  if (length(x) != length(order)) stop("length of order has to equal the length of x")

  out <- setNames(x, paste0("beta", order))
  returnPWMs(out, order, input = x, func = "as.PWMs")
}

#' @describeIn as.PWMs as.PWMs for numeric data matrices
#' @method as.PWMs matrix
#' @export
as.PWMs.matrix <- function(x, order = row(x)[, 1]-1, ...) {
  if (nrow(x) != length(order)) stop("length of order has to equal the number of rows of x")

  out <- x
  rownames(out) <- paste0("beta", order)
  returnPWMs(out, order, input = x, func = "as.PWMs")
}

#' @describeIn as.PWMs as.PWMs for numeric data lists
#' @method as.PWMs list
#' @export
as.PWMs.list <- function(x, order = seq_along(x[[1]])-1, ...) {
  if (any(vapply(x, length, numeric(1)) != length(order))) stop("all elements have to be of the length of order. ")

  out <- x
  for (i in 1L:length(out)) names(out[[i]]) <- paste0("beta", order)

  returnPWMs(out, order, input = x, func = "as.PWMs")
}

#' @describeIn as.PWMs as.PWMs for numeric data.frames
#' @method as.PWMs data.frame
#' @export
as.PWMs.data.frame <- function(x, formula, order = NULL, ...) {

  nam <- getFormulaSides(formula, names(x))

  if (is.null(order)) order <- seq_along(nam$lhs)-1
  if (length(nam$lhs) != length(order)) stop("all elements have to be of the length of order. ")

  out <- cbind(
    x[nam$rhs],
    as.data.frame(t(apply(x[nam$lhs], 1, as.PWMs, order = order)))
  )

  returnPWMs(out, order, input = x, formula = formula, func = "as.PWMs")
}
