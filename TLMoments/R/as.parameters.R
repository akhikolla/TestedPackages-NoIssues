#' @title
#' Converting to parameters-objects
#' @description
#' Convert vector, matrix, list, or data.frame to parameters-objects.
#'
#' @param ... parameters of distribution. This can be named vectors or lists, matrices, or data.frames. See examples below.
#' @param x numeric vector, matrix, list, or data.frame of parameters.
#' @param formula if \code{x} is data.frame a formula has to be given.
#' @param distr character giving the distribution. Function of name
#' q\"distr\" has to be available.
#'
#' @return object of class parameters, see parameters help page.
#' @seealso \code{\link{parameters}}
#' @examples
#' # Vector input:
#' as.parameters(loc = 3, scale = 2, shape = .4, distr = "gev")
#' as.parameters(c(loc = 3, scale = 2, shape = .4), distr = "gev")
#'
#' # Names can be shortened if unambiguous:
#' as.parameters(l = 3, sc = 2, sh = .4, distr = "gev")
#' as.parameters(m = 3, s = 1, distr = "norm")
#'
#' # Wrong or ambiguous names lead to errors!
#' \dontrun{
#' as.parameters(l = 3, s = 2, s = .4, distr = "gev")
#' as.parameters(loc2 = 3, scale = 2, shape = .4, distr = "gev")
#' }
#'
#' # If no names are given, a warning is given and they are guessed for gev, gpd, gum, and ln3.
#' as.parameters(3, 2, .4, distr = "gev")
#' as.parameters(c(3, 2, .4), distr = "gev")
#' \dontrun{
#' as.parameters(3, 2, .2, .4, distr = "gev") #=> doesn't work
#' }
#'
#' # Matrix input:
#' # Parameters in matrices must have either matching rownames or colnames!
#' as.parameters(cbind(loc = 10, scale = 4, shape = seq(0, .4, .1)), distr = "gev")
#' as.parameters(rbind(loc = 10, scale = 4, shape = seq(0, .4, .1)), distr = "ln3")
#'
#' # If no names are given, a guess is made based on number of rows
#' # or cols according to distribution (and a warning is given).
#' as.parameters(matrix(1:9, nr = 3), distr = "gev")
#' as.parameters(matrix(1:8, nc = 2), distr = "gum")
#'
#'
#' # The same principles apply for list input and data.frames:
#'
#' # List input:
#' as.parameters(list(list(mean = 2, sd = 1), list(mean = 0, sd = 1)), distr = "norm")
#' as.parameters(list(c(m = 2, s = 1), c(m = 0, s = 1)), distr = "norm")
#' as.parameters(list(c(loc = 2, scale = 1), c(0, 1)), distr = "gum")
#' \dontrun{
#' as.parameters(list(c(loc = 2, scale = 1), c(0, 1, 2)), distr = "gum")
#' }
#'
#' # Dataframe input:
#' xdat <- data.frame(station = c(1, 2), mean = c(2, 0), sd = c(1, 1))
#' as.parameters(xdat, cbind(mean, sd) ~ station, distr = "norm")
#' as.parameters(xdat, . ~ station, distr = "norm")
#' as.parameters(xdat, cbind(mean, sd) ~ ., distr = "norm")
#'
#' xdat <- data.frame(station = c(1, 2), m = c(2, 0), s = c(1, 1))
#' as.parameters(xdat, cbind(m, s) ~ station, distr = "norm")
#' \dontrun{
#' as.parameters(xdat, cbind(m, s) ~ station, distr = "gev")
#' }
#'
#' ###
#'
#' # Results of as.parameters can be used in the normal TLMoments-scheme:
#' # they can be transfered to quantiles or to TLMoments.
#'
#' xdat <- data.frame(station = c(1, 2), mean = c(2, 0), sd = c(1, 1))
#' quantiles(as.parameters(xdat, cbind(mean, sd) ~ ., distr = "norm"), c(.99))
#'
#' # quantile estimation
#' p <- as.parameters(loc = 3, scale = 2, shape = .4, distr = "gev")
#' quantiles(p, c(.9, .95))
#' p <- as.parameters(cbind(loc = 10, scale = 4, shape = seq(0, .4, .1)), distr = "gev")
#' quantiles(p, c(.9, .95))
#' p <- as.parameters(list(list(mean = 2, sd = 1), list(mean = 0, sd = 1)), distr = "norm")
#' quantiles(p, c(.95, .975))
#'
#' # With magrittr
#' library(magrittr)
#' as.parameters(loc = 3, scale = 2, shape = .4, distr = "gev") %>% quantiles(c(.9, .99))
#' @rdname as.parameters
#' @export
as.parameters <- function(..., distr = NULL) {
  if (is.null(distr)) stop("distr of distribution has to be submitted")

  # if (grepl("::", x = distr)) { # Falls pkg::func
  #   f <- sub("^([a-zA-Z0-9]*)::([a-zA-Z0-9]*)$", "\\1::q\\2", x = distr)
  #   q <- eval(parse(text = paste0("match.fun(", f,")")))
  # } else { # falls nur func
  #   q <- eval(parse(text = paste0("match.fun(q", distr, ")")))
  # }
  # if (!is.function(q)) stop(paste0("Found no q-function for ", distr))

  UseMethod("as.parameters")
}

#' @describeIn as.parameters as.parameters for numeric data vectors
#' @method as.parameters numeric
#' @export
as.parameters.numeric <- function(..., distr) {
  out <- unlist(list(...))

  # Assume names if not available. Only for GEV, GPD, GUM, LN3
  if (is.null(names(out))) {
    if (distr %in% c("gev", "gpd", "ln3") & length(out) == 3) {
      names(out) <- c("loc", "scale", "shape")
      warning("Non-named parameters. Assuming loc, scale, and shape. ")
    } else if (distr %in% c("gum") & length(out) == 2) {
      names(out) <- c("loc", "scale")
      warning("Non-named parameters. Assuming loc and scale. ")
    } else {
      stop("Please provide parameter names. ")
    }
  }

  names(out) <- checkParameterNames(names(out), distr)

  returnParameters(
    out, distr,
    func = "as.parameters",
    input = out
  )
}

#' @describeIn as.parameters as.parameters for numeric data matrices
#' @method as.parameters matrix
#' @export
as.parameters.matrix <- function(x, distr, ...) {

  if (is.null(dimnames(x))) dimnames(x) <- list(NULL, NULL)
  dimns <- lapply(dimnames(x), function(xx) {
    tryCatch(checkParameterNames(xx, distr), error = function(e) NULL)
  })

  .d <- sapply(dimns, is.null)

  # Both colnames and rownames
  if (all(.d == c(FALSE, FALSE))) {
    stop("Both rownames and colnames contain parameter names. ")
  }

  # Only colnames
  if (all(.d == c(TRUE, FALSE))) {
    x <- t(x)
  }

  # No names at all:
  if (all(.d == c(TRUE, TRUE))) {
    if (distr %in% c("gev", "gpd", "ln3")) {
      if (nrow(x) == 3) {
        rownames(x) <- c("loc", "scale", "shape")
        warning("Non-named parameters. Assuming loc, scale, and shape. ")
      } else if (ncol(x) == 3) {
        x <- t(x)
        rownames(x) <- c("loc", "scale", "shape")
        warning("Non-named parameters. Assuming loc, scale, and shape. ")
      } else {
        stop("Please provide parameter names. ")
      }

    } else if (distr == "gum") {
      if (nrow(x) == 2) {
        rownames(x) <- c("loc", "scale")
        warning("Non-named parameters. Assuming loc and scale. ")
      } else if (ncol(x) == 2) {
        x <- t(x)
        rownames(x) <- c("loc", "scale")
        warning("Non-named parameters. Assuming loc and scale. ")
      } else {
        stop("Please provide parameter names. ")
      }
    }
  }

  returnParameters(
    x, distr,
    func = "as.parameters",
    input = x
  )
}

#' @describeIn as.parameters as.parameters for numeric data lists
#' @method as.parameters list
#' @export
as.parameters.list <- function(x, distr, ...) {
  out <- lapply(x, function(y) {
    do.call(as.parameters, args = list(unlist(y), distr = distr))
  })

  # Delete attributes...
  for (i in seq_along(out)) {
    attr(out[[i]], "source") <- NULL
    attr(out[[i]], "distribution") <- NULL
    class(out[[i]]) <- class(out[[i]])[-1]
  }
  # ...and add global attributes
  returnParameters(
    out, distr,
    func = "as.parameters",
    input = out
  )
}

#' @describeIn as.parameters as.parameters for numeric data.frames
#' @method as.parameters data.frame
#' @export
as.parameters.data.frame <- function(x, formula, distr, ...) {

  nam <- getFormulaSides(formula, names(x))
  out <- cbind(x[nam$rhs], x[nam$lhs])
  names(out)[names(out) %in% nam$lhs] <- checkParameterNames(nam$lhs, distr)

  returnParameters(
    out, distr,
    func = "as.parameters",
    input = x,
    formula = nam$new.formula
  )
}
