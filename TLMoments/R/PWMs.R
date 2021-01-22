#' @title
#' Probability weighted moments
#' @description
#' Calculates probability weighted moments up to a specific order. Note that PWMs start with
#' order 0. Acceptable input types are numeric vectors, matrices, lists, and data.frames.
#'
#' @param x numeric vector or matrix, list, or data.frame of data OR an object of TLMoments.
#' @param formula if x is of type data.frame a formula has to be submitted.
#' @param max.order integer, maximal order of PWMs.
#' @param na.rm logical, indicates if NAs should be removed.
#' @param ... additional arguments.
#'
#' @return numeric vector, matrix, list, or data.frame consisting of the PWMs and
#' with class \code{PWMs}.
#' The object contains the following attributes: \itemize{
#'  \item \code{order}: a integer vector with corresponding PWM orders
#'  \item \code{source}: a list with background information (used function, data, n, formula;
#'  mainly for internal purposes)
#' }
#' The attributes are hidden in the print-function for a clearer presentation.
#' @references Greenwood, J. A., Landwehr, J. M., Matalas, N. C., & Wallis, J. R. (1979). Probability weighted moments: definition and relation to parameters of several distributions expressable in inverse form. Water Resources Research, 15(5), 1049-1054.
#'
#' @examples
#' # Generating data sets:
#' xmat <- matrix(rnorm(100), nc = 4)
#' xvec <- xmat[, 3]
#' xlist <- lapply(1L:ncol(xmat), function(i) xmat[, i])
#' xdat <- data.frame(
#'  station = rep(letters[1:2], each = 50),
#'  season = rep(c("S", "W"), 50),
#'  hq = as.vector(xmat)
#' )
#'
#' # Calculating PWMs from data:
#' PWMs(xvec)
#' PWMs(xmat)
#' PWMs(xlist)
#' PWMs(xdat, formula = hq ~ station)
#' PWMs(xdat, formula = hq ~ season)
#' PWMs(xdat, formula = hq ~ .)
#' PWMs(xdat, formula = . ~ station + season)
#'
#' # Calculating PWMs from L-moments:
#' PWMs(TLMoments(xvec))
#' PWMs(TLMoments(xmat))
#' PWMs(TLMoments(xlist))
#' PWMs(TLMoments(xdat, hq ~ station))
#' PWMs(TLMoments(xdat, hq ~ season))
#' PWMs(TLMoments(xdat, hq ~ .))
#' PWMs(TLMoments(xdat, . ~ station + season))
#'
#' # In data.frame-mode invalid names are preceded by "."
#' xdat <- data.frame(
#'  beta0 = rep(letters[1:2], each = 50),
#'  beta1 = as.vector(xmat)
#' )
#' PWMs(xdat, formula = beta1 ~ beta0)
#'
#' @rdname PWMs
#' @export
PWMs <- function(x, ...) UseMethod("PWMs")


#' @title returnPWMs
#' @description Sets attributes to PWMs objects and returns them. This function is for internal use.
#' @param out -
#' @param order -
#' @param ... -
#' @return An object of class PWMs.
returnPWMs <- function(out, order, ...) {

  class <- class(out)
  args <- list(...)

  # If no func attribute is set, set to
  if (!exists("func", args)) args$func <- "PWMs"

  # If more than one func attributes are given, concatenate them
  if (sum(names(args) == "func") >= 2) {
    newfunc <- as.vector(unlist(args[names(args) == "func"]))
    args$func <- NULL
    args$func <- newfunc
  }

  # Calculate n if not available and data exists.
  if (!exists("n", args) && exists("data", args)) {# && args$func == "PWMs") {
    if (inherits(out, "numeric")) {
      args$n <- sum(!is.na(args$data))
    } else if (inherits(out, "matrix")) {
      args$n <- apply(args$data, 2, function(y) sum(!is.na(y)))
    } else if (inherits(out, "list")) {
      args$n <- vapply(args$data, length, numeric(1))
    } else if (inherits(out, "data.frame")) {
      args$n <- aggregate(args$formula, args$data, length)[[getFormulaSides(args$formula)$lhs]]
    }
  }

  # Attributes of PWMs
  # order
  # source: func
  #         data (if calculated)
  #         n (if calculated)
  #         formula (if data is data.frame)
  #         lambdas (if coming from TLMoments)
  #         trimmings (if coming from TLMoments)
  # class: "PWMs"

  attr(out, "order") <- order
  attr(out, "source") <- args
  class(out) <- c("PWMs", class)

  out
}

#' @rdname PWMs
#' @method PWMs numeric
#' @export
PWMs.numeric <- function(x, max.order = 4L, na.rm = FALSE, ...) {
  out <- setNames(
    PWM(x, order = 0L:max.order, na.rm = na.rm),
    paste0("beta", 0:max.order)
  )

  returnPWMs(out, 0L:max.order, data = x)
}

#' @rdname PWMs
#' @method PWMs matrix
#' @export
PWMs.matrix <- function(x, max.order = 4L, na.rm = FALSE, ...) {
  out <- apply(x, 2, PWM, order = 0L:max.order, na.rm = na.rm)

  returnPWMs(out, 0L:max.order, data = x)
}

#' @rdname PWMs
#' @method PWMs list
#' @export
PWMs.list <- function(x, max.order = 4L, na.rm = FALSE, ...) {
  out <- lapply(x, PWM, order = 0L:max.order, na.rm = na.rm)

  returnPWMs(out, 0L:max.order, data = x)
}

#' @rdname PWMs
#' @method PWMs data.frame
#' @export
PWMs.data.frame <- function(x, formula, max.order = 4L, na.rm = FALSE, ...) {

  # Check for and repair invalid variables named beta[0-9]
  x <- correctNames(x, "beta[0-9]*", ".")
  formula <- correctNames(formula, "beta[0-9]*", ".")

  nam <- getFormulaSides(formula, names(x))
  r <- aggregate(nam$new.formula, data = x, FUN = PWM, order = 0L:max.order, na.rm = na.rm)
  out <- cbind(r[-length(r)], as.data.frame(r[[length(r)]]))

  returnPWMs(out, 0L:max.order, data = x, formula = nam$new.formula)
}


#' @rdname PWMs
#' @method PWMs TLMoments
#' @export
PWMs.TLMoments <- function(x, ...) {
  if (attr(x, "leftrim") != 0 | attr(x, "rightrim") != 0) stop("Transformation to PWMs only works for L-moments. ")
  if (any(diff(attr(x, "order"))) != 1) stop("Transformation to PWMs only runs for sequent L-moments. ")

   UseMethod("PWMs.TLMoments", x$lambdas)
}

#' @method PWMs.TLMoments numeric
#' @export
PWMs.TLMoments.numeric <- function(x, ...) {
  max.order <- max(attr(x, "order"))
  out <- as.numeric(solve(Z_C(max.order, 0, 0)) %*% x$lambdas)
  names(out) <- paste0("beta", 0:(max.order-1))

  do.call(
    returnPWMs, c(
      list(out = out, order = 0L:(max.order-1)),
      func = "PWMs",
      lambdas = list(x$lambdas),
      trimmings = list(c(attr(x, "leftrim"), attr(x, "rightrim"))),
      tl.order = list(attr(x, "order")),
      attr(x, "source")
    )
  )
}

#' @method PWMs.TLMoments matrix
#' @export
PWMs.TLMoments.matrix <- function(x, ...) {
  max.order <- max(attr(x, "order"))
  out <- solve(Z_C(max.order, 0, 0)) %*% x$lambdas
  rownames(out) <- paste0("beta", 0:(max.order-1))

  do.call(
    returnPWMs, c(
      list(out = out, order = 0L:(max.order-1)),
      func = "PWMs",
      lambdas = list(x$lambdas),
      trimmings = list(c(attr(x, "leftrim"), attr(x, "rightrim"))),
      tl.order = list(attr(x, "order")),
      attr(x, "source")
    )
  )
}

#' @method PWMs.TLMoments list
#' @export
PWMs.TLMoments.list <- function(x, ...) {
  max.order <- max(attr(x, "order"))
  A <- solve(Z_C(max.order, 0, 0))
  out <- lapply(x$lambdas, function(x) setNames(as.numeric(A %*% x), paste0("beta", 0:(max.order-1))))

  do.call(
    returnPWMs, c(
      list(out = out, order = 0L:(max.order-1)),
      func = "PWMs",
      lambdas = list(x$lambdas),
      trimmings = list(c(attr(x, "leftrim"), attr(x, "rightrim"))),
      tl.order = list(attr(x, "order")),
      attr(x, "source")
    )
  )
}

#' @method PWMs.TLMoments data.frame
#' @export
PWMs.TLMoments.data.frame <- function(x, ...) {
  max.order <- max(attr(x, "order"))

  ls <- x$lambdas[, grep("L[0-9]*", names(x$lambdas))]
  fac <- x$lambdas[, !grepl("L[0-9]*", names(x$lambdas)), drop = FALSE]

  out <- as.data.frame(
    t(solve(Z_C(max.order, 0, 0)) %*% t(as.matrix(ls)))
  )
  names(out) <- paste0("beta", 0:(max.order-1))
  out <- cbind(fac, out)

  do.call(
    returnPWMs, c(
      list(out = out, order = 0L:(max.order-1)),
      func = "PWMs",
      lambdas = list(x$lambdas),
      trimmings = list(c(attr(x, "leftrim"), attr(x, "rightrim"))),
      tl.order = list(attr(x, "order")),
      attr(x, "source")
    )
  )
}

#' @export
print.PWMs <- function(x, ...) {
  if (inherits(x, "data.frame")) {
    print.data.frame(x)
    return(invisible(x))
  }

  tmp <- x
  attributes(tmp) <- NULL

  dim(tmp) <- dim(x)
  names(tmp) <- names(x)
  dimnames(tmp) <- dimnames(x)

  print(tmp)
  invisible(x)
}
