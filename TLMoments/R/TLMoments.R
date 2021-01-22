#' @title
#' Trimmed L-moments
#' @description
#' Calculates empirical or theoretical Trimmed L-moments and -ratios up to a specific order.
#' If empirical moments should be calculated, acceptable input types are numeric vectors,
#' matrices, lists, data.frames. TLMoments is type-preservative, so the input type is also
#' the output type. If theoretical moments should be calculated, the input type has to be
#' of class parameters or PWMs, so an object returned by parameters, as.parameters or
#' PWMs, as.PWMs.
#' @param x numeric data in form of vector, matrix, list, or data.frame OR an object
#' of class parameters or PWMs.
#' @param formula if \code{x} is data.frame. See examples.
#' @param leftrim integer indicating lower trimming parameter, has to be greater than 0.
#' @param rightrim integer indicating upper trimming parameter, has to be greater than 0.
#' @param max.order integer, maximum order of Trimmed L-moments/ratios, has to be
#' greater than 1.
#' @param na.rm logical, indicates if NAs should be removed. Only if empirical moments
#' are calculated.
#' @param computation.method character, indicating if the computation is performed via
#' PWMs, direct, recursive, or recurrence (see References Hosking & Balakrishnan, 2015).
#' Possible values are \code{auto} (default, automatically choose appropriate method), \code{pwm},
#' \code{direct}, \code{recursive}, or \code{recurrence}. Only if empirical moments are calculated.
#' @param ... additional arguments.
#'
#' @return list of two dimensions: \code{lambdas}/\code{ratios} are a numeric vector, matrix,
#' list, or data.frame consisting of the TL-moments/TL-moment-ratios. The list has the class
#' \code{TLMoments}.
#' The object contains the following attributes: \itemize{
#'  \item \code{leftrim}: a numeric giving the used leftrim-argument
#'  \item \code{rightrim}: a numeric giving the used rightrim-argument
#'  \item \code{order}: a integer vector with corresponding TL-moment orders
#'  \item \code{source}: a list with background information (used function, data, n, formula, computation method;
#'  mainly for internal purposes)
#' }
#' The attributes are hidden in the print-function for a clearer presentation.
#'
#' @references Elamir, E. A., & Seheult, A. H. (2003). Trimmed L-moments. Computational Statistics & Data Analysis, 43(3), 299-314.
#' @references Hosking, J. R. M. (1990). L-moments: analysis and estimation of distributions using linear combinations of order statistics. Journal of the Royal Statistical Society. Series B (Methodological), 105-124.
#' @references Hosking, J. R. M. (2007). Some theory and practical uses of trimmed L-moments. Journal of Statistical Planning and Inference, 137(9), 3024-3039.
#' @references Hosking, J. R. M., & Balakrishnan, N. (2015). A uniqueness result for L-estimators, with applications to L-moments. Statistical Methodology, 24, 69-80.
#' @seealso \code{\link{PWMs}}, \code{\link{parameters}}, \code{\link{quantiles}}, \code{\link{summary.TLMoments}}, \code{\link{as.TLMoments}}
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
#' # Calculating TL-moments from data:
#' TLMoments(xvec, leftrim = 0, rightrim = 1)
#' TLMoments(xmat, leftrim = 1, rightrim = 1)
#' TLMoments(xlist, max.order = 7)
#' TLMoments(xdat, hq ~ station, leftrim = 0, rightrim = 2)
#' TLMoments(xdat, hq ~ season, leftrim = 0, rightrim = 2)
#' TLMoments(xdat, hq ~ ., leftrim = 0, rightrim = 2)
#'
#' # Calculating TL-moments from PWMs:
#' TLMoments(PWMs(xvec))
#' TLMoments(PWMs(xmat), rightrim = 1)
#' TLMoments(PWMs(xlist), leftrim = 1, rightrim = 1)
#' TLMoments(PWMs(xdat, hq ~ station), leftrim = 0, rightrim = 2)
#' TLMoments(PWMs(xdat, hq ~ station + season), leftrim = 0, rightrim = 2)
#' TLMoments(as.PWMs(cbind(c(0.12, .41, .38, .33), c(.05, 0.28, .25, .22))), 0, 1)
#'
#' # Calculating TL-moments from parameters:
#' (tlm <- TLMoments(xmat, leftrim = 0, rightrim = 1))
#' TLMoments(parameters(tlm, "gev"))
#'
#' (tlm <- TLMoments(xdat, hq ~ station, leftrim = 0, rightrim = 2))
#' TLMoments(parameters(tlm, "gev"))
#'
#' p <- as.parameters(loc = 3, scale = 2, shape = .4, distr = "gev")
#' TLMoments(p, rightrim = 1)
#'
#' p <- as.parameters(cbind(loc = 10, scale = 4, shape = seq(0, .4, .1)), distr = "gev")
#' TLMoments(p, max.order = 6)
#'
#' p <- as.parameters(list(
#'  list(loc = 3, scale = 2, shape = .4),
#'  list(loc = 3, scale = 2, shape = .2)
#' ), distr = "gev")
#' TLMoments(p)
#'
#' p <- as.parameters(data.frame(
#'  station = letters[1:2],
#'  loc = c(2, 3),
#'  scale = c(2, 2),
#'  shape = c(.4, .2)
#' ), .~station, distr = "gev")
#' TLMoments(p)
#' @export
TLMoments <- function(x, ...) {

  args <- list(...)
  # if (exists("leftrim", args) && exists("rightrim", args))
  #  if (!are.integer.like(args$leftrim, args$rightrim) | any(c(args$leftrim, args$rightrim) < 0))
  #   stop("leftrim and rightrim must be positive integers. ")
  if ("max.order" %in% names(args) && !are.integer.like(args$max.order))
    stop("max.order must be integer-like. ")
  if ("na.rm" %in% names(args) && !is.logical(args$na.rm))
    stop("na.rm must be TRUE or FALSE. ")

  UseMethod("TLMoments")
}


#' @title returnTLMoments
#' @description Sets attributions to TLMoments objects and returns them. This function is for internal use.
#' @param out -
#' @param leftrim -
#' @param rightrim -
#' @param order -
#' @param ... -
#' @return An object of class TLMoments.
returnTLMoments <- function(out, leftrim, rightrim, order, ...) {

  class <- class(out$lambdas)
  args <- list(...)

  # If no func attribute is set, set to
  if (!exists("func", args)) args$func <- "TLMoments"

  # If more than one func attributes are given, concatenate them
  if (sum(names(args) == "func") >= 2) {
    newfunc <- as.vector(unlist(args[names(args) == "func"]))
    args$func <- NULL
    args$func <- newfunc
  }

  # Calculate n if not available and data exists.
  if (!exists("n", args) && exists("data", args)) {
    if (inherits(out$lambdas, "numeric")) {
      args$n <- sum(!is.na(args$data))
    } else if (inherits(out$lambdas, "matrix")) {
      args$n <- apply(args$data, 2, function(y) sum(!is.na(y)))
    } else if (inherits(out$lambdas, "list")) {
      args$n <- vapply(args$data, length, numeric(1))
    } else if (inherits(out$lambdas, "data.frame")) {
      args$n <- aggregate(args$formula, args$data, length)[[getFormulaSides(args$formula)$lhs]]
    }
  }

  # Attributes of TLMoments
  # leftrim
  # rightrim
  # order
  # source: func
  #         computation.method (if calculated)
  #         data (if calculated)
  #         input (if not calculated)
  #         n (if calculated)
  #         formula (if data is data.frame)
  #         pwms (if coming from PWMs)
  # class: "TLMoments", "list"

  attr(out, "leftrim") <- leftrim
  attr(out, "rightrim") <- rightrim
  attr(out, "order") <- order
  attr(out, "source") <- args
  class(out) <- c("TLMoments", "list")

  out
}


#' @describeIn TLMoments TLMoments for numeric vector of data
#' @method TLMoments numeric
#' @export
TLMoments.numeric <- function(x, leftrim = 0L, rightrim = 0L, max.order = 4L,
                              na.rm = FALSE, computation.method = "auto",
                              ...) {
  if (computation.method == "auto")
    computation.method <- select_computation(leftrim, rightrim)

  ls <- TLMoment(x, order = 1L:max.order, leftrim = leftrim, rightrim = rightrim, na.rm = na.rm, computation.method = computation.method)
  out <- list(
    lambdas = ls,
    ratios = calcRatios(ls)
  )

  returnTLMoments(out, leftrim, rightrim, 1L:max.order,
                  tl.computation.method = computation.method, data = x)
}

#' @describeIn TLMoments TLMoments for numeric matrix of data
#' @method TLMoments matrix
#' @export
TLMoments.matrix <- function(x, leftrim = 0L, rightrim = 0L, max.order = 4L,
                             na.rm = FALSE, computation.method = "auto",
                             ...) {
  if (computation.method == "auto")
    computation.method <- select_computation(leftrim, rightrim)

  ls <- apply(x, 2, TLMoment, order = 1L:max.order, leftrim = leftrim, rightrim = rightrim, na.rm = na.rm, computation.method = computation.method)
  if (max.order == 1) {
    dim(ls) <- c(1, ncol(x))
  }
  out <- list(
    lambdas = ls,
    ratios = apply(ls, 2, calcRatios)
  )

  returnTLMoments(out, leftrim, rightrim, 1L:max.order,
                  tl.computation.method = computation.method, data = x)
}

#' @describeIn TLMoments TLMoments for numeric list of data
#' @method TLMoments list
#' @export
TLMoments.list <- function(x, leftrim = 0L, rightrim = 0L, max.order = 4L, na.rm = FALSE,
                           computation.method = "auto",
                           ...) {
  if (computation.method == "auto")
    computation.method <- select_computation(leftrim, rightrim)

  ls <- lapply(x, TLMoment, order = 1L:max.order, leftrim = leftrim, rightrim = rightrim, na.rm = na.rm, computation.method = computation.method)
  out <- list(
    lambdas = ls,
    ratios = lapply(ls, calcRatios)
  )

  returnTLMoments(out, leftrim, rightrim, 1L:max.order,
                  tl.computation.method = computation.method, data = x)
}

#' @describeIn TLMoments TLMoments for numeric data.frame of data
#' @method TLMoments data.frame
#' @export
TLMoments.data.frame <- function(x, formula, leftrim = 0L, rightrim = 0L, max.order = 4L,
                                 na.rm = FALSE, computation.method = "auto",
                                 ...) {
  if (computation.method == "auto")
    computation.method <- select_computation(leftrim, rightrim)

  # Check for and repair invalid variables named [L|T][0-9]
  x <- correctNames(x, "[L|T][0-9]*", ".")
  formula <- correctNames(formula, "[L|T][0-9]*", ".")

  nam <- getFormulaSides(formula, names(x))
  agg <- aggregate(nam$new.formula, data = x, FUN = TLMoment, order = 1L:max.order, leftrim = leftrim, rightrim = rightrim, na.rm = na.rm, computation.method = computation.method)
  out <- list(
    lambdas = cbind(agg[-length(agg)], as.data.frame(agg[[length(agg)]])),
    ratios = cbind(agg[-length(agg)], as.data.frame(t(apply(agg[[length(agg)]], 1, calcRatios)[-1, ])))
  )

  returnTLMoments(out, leftrim, rightrim, 1L:max.order,
                  tl.computation.method = computation.method, data = x, formula = nam$new.formula)
}



#' @describeIn TLMoments TLMoments for PWMs-object
#' @method TLMoments PWMs
#' @export
TLMoments.PWMs <- function(x, leftrim = 0L, rightrim = 0L, ...) {
  if (any(diff(attr(x, "order")) != 1))
    stop("PWM order must not have gaps")

  UseMethod("TLMoments.PWMs")
}

#' @method TLMoments.PWMs numeric
#' @export
TLMoments.PWMs.numeric <- function(x, leftrim = 0L, rightrim = 0L, ...) {
  ls <- PWM_to_TLMoments(x, leftrim, rightrim)
  out <- list(
    lambdas = setNames(ls, paste0("L", seq_along(ls))),
    ratios = calcRatios(ls)
  )

  do.call(returnTLMoments, c(
    list(out = out, leftrim = leftrim, rightrim = rightrim, order = seq_along(ls)),
    func = "TLMoments.PWMs",
    pwms = list(removeAttributes(x)),
    attr(x, "source")
  ))
}

#' @method TLMoments.PWMs matrix
#' @export
TLMoments.PWMs.matrix <- function(x, leftrim = 0L, rightrim = 0L, ...) {
  ls <- apply(x, 2, function(xx) {
    PWM_to_TLMoments(xx, leftrim, rightrim)
  })
  rownames(ls) <- paste0("L", seq_along(ls[, 1]))
  out <- list(
    lambdas = ls,
    ratios = apply(ls, 2, calcRatios)
  )

  do.call(returnTLMoments, c(
    list(out = out, leftrim = leftrim, rightrim = rightrim, order = seq_along(ls[, 1])),
    func = "TLMoments.PWMs",
    pwms = list(removeAttributes(x)),
    attr(x, "source")
  ))
}

#' @method TLMoments.PWMs list
#' @export
TLMoments.PWMs.list <- function(x, leftrim = 0L, rightrim = 0L, ...) {
  ls <- lapply(x, function(xx) {
    PWM_to_TLMoments(xx, leftrim, rightrim)
  })
  ls <- lapply(ls, function(l) setNames(l, paste0("L", seq_along(l))))
  out <- list(
    lambdas = ls,
    ratios = lapply(ls, calcRatios)
  )

  do.call(returnTLMoments, c(
    list(out = out, leftrim = leftrim, rightrim = rightrim, order = seq_along(ls[[1]])),
    func = "TLMoments.PWMs",
    pwms = list(removeAttributes(x)),
    attr(x, "source")
  ))
}

#' @method TLMoments.PWMs data.frame
#' @export
TLMoments.PWMs.data.frame <- function(x, leftrim = 0L, rightrim = 0L, ...) {
  pwms <- x[, grep("beta[0-9]*", names(x)), drop = FALSE]
  fac <- x[, !grepl("beta[0-9]*", names(x)), drop = FALSE]
  ls <- apply(pwms, 1, function(xx) {
    PWM_to_TLMoments(xx, leftrim, rightrim)
  })
  ratios <- apply(ls, 2, calcRatios)
  ratios <- as.data.frame(t(ratios))
  lambdas <- as.data.frame(t(ls))
  names(lambdas) <- paste0("L", 1L:ncol(lambdas))

  out <- list(
    lambdas = cbind(fac, lambdas),
    ratios = cbind(fac, ratios)
  )

  do.call(returnTLMoments, c(
    list(out = out, leftrim = leftrim, rightrim = rightrim, order = 1L:ncol(lambdas)),
    func = "TLMoments.PWMs",
    pwms = list(removeAttributes(x)),
    attr(x, "source")
  ))
}



#' @describeIn TLMoments TLMoments for parameters-object
#' @method TLMoments parameters
#' @export
TLMoments.parameters <- function(x,
                                 leftrim = attr(x, "source")$trimmings[1],
                                 rightrim = attr(x, "source")$trimmings[2],
                                 max.order = 4L,
                                 ...) {

  if (!is.null(max.order) && !are.integer.like(max.order))
    stop("max.order must be integer-like. ")
  if (!is.null(max.order) && max.order <= 0)
    stop("max.order must be positive. ")
  if ((!is.null(leftrim) && !are.integer.like(leftrim)) || (!is.null(rightrim) && !are.integer.like(rightrim)))
    stop("leftrim and rightrim must be integer-like. ")
  if ((!is.null(leftrim) && leftrim < 0) || (!is.null(rightrim) && rightrim < 0))
    stop("leftrim and rightrim must be non-negative. ")

  UseMethod("TLMoments.parameters")
}

#' @method TLMoments.parameters numeric
#' @export
TLMoments.parameters.numeric <- function(x,
                                         leftrim = attr(x, "source")$trimmings[1],
                                         rightrim = attr(x, "source")$trimmings[2],
                                         max.order = 4L,
                                         ...) {

  if (is.null(leftrim)) leftrim <- 0L
  if (is.null(rightrim)) rightrim <- 0L

  if (is.null(attr(x, "source")$lambdas) ||
      !identical(attr(x, "source")$trimmings, c(leftrim, rightrim)) ||
      (!is.null(attr(x, "source")$tl.order) && max(attr(x, "source")$tl.order) != max.order)) { # calculate new

    ls <- calcTLMom(max.order, leftrim, rightrim, qfunc = getQ(x))
  } else { # or use old calculations
    ls <- attr(x, "source")$lambdas
  }

  out <- list(
    lambdas = setNames(ls, paste0("L", seq_along(ls))),
    ratios = calcRatios(ls)
  )

  do.call(returnTLMoments,c(
    list(out = out, leftrim = leftrim, rightrim = rightrim, order = 1L:max.order),
    func = "TLMoments.parameters",
    distr = attr(x, "distribution"),
    parameters = list(removeAttributes(x)),
    attr(x, "source")
  ))
}

#' @method TLMoments.parameters matrix
#' @export
TLMoments.parameters.matrix <- function(x,
                                        leftrim = attr(x, "source")$trimmings[1],
                                        rightrim = attr(x, "source")$trimmings[2],
                                        max.order = 4L,
                                        ...) {

  if (is.null(leftrim)) leftrim <- 0L
  if (is.null(rightrim)) rightrim <- 0L

  if (is.null(attr(x, "source")$lambdas) ||
      !identical(attr(x, "source")$trimmings, c(leftrim, rightrim)) ||
      (!is.null(attr(x, "source")$tl.order) && max(attr(x, "source")$tl.order) != max.order)) { # calculate new

    ls <- apply(x, 2, function(xx) {
      calcTLMom(max.order, leftrim, rightrim,
                qfunc = do.call(getQ, c(x = attr(x, "distribution"), as.list(xx))))
    })
  } else { # or use old calculations
    ls <- attr(x, "source")$lambdas
  }

  rownames(ls) <- paste0("L", seq_along(ls[, 1]))
  out <- list(
    lambdas = ls,
    ratios = apply(ls, 2, calcRatios)
  )

  do.call(returnTLMoments,c(
    list(out = out, leftrim = leftrim, rightrim = rightrim, order = 1L:max.order),
    func = "TLMoments.parameters",
    distr = attr(x, "distribution"),
    parameters = list(removeAttributes(x)),
    attr(x, "source")
  ))
}

#' @method TLMoments.parameters list
#' @export
TLMoments.parameters.list <- function(x,
                                      leftrim = attr(x, "source")$trimmings[1],
                                      rightrim = attr(x, "source")$trimmings[2],
                                      max.order = 4L,
                                      ...) {

  if (is.null(leftrim)) leftrim <- 0L
  if (is.null(rightrim)) rightrim <- 0L

  if (is.null(attr(x, "source")$lambdas) ||
      !identical(attr(x, "source")$trimmings, c(leftrim, rightrim)) ||
      (!is.null(attr(x, "source")$tl.order) && max(attr(x, "source")$tl.order) != max.order)) { # calculate new

    ls <- lapply(x, function(xx) {
      calcTLMom(max.order, leftrim, rightrim,
                qfunc = do.call(getQ, c(x = attr(x, "distribution"), as.list(xx))))
    })
  } else { # or use old calculations
    ls <- attr(x, "source")$lambdas
  }

  ls <- lapply(ls, function(l) setNames(l, paste0("L", seq_along(l))))
  out <- list(
    lambdas = ls,
    ratios = lapply(ls, calcRatios)
  )

  do.call(returnTLMoments,c(
    list(out = out, leftrim = leftrim, rightrim = rightrim, order = 1L:max.order),
    func = "TLMoments.parameters",
    distr = attr(x, "distribution"),
    parameters = list(removeAttributes(x)),
    attr(x, "source")
  ))
}

#' @method TLMoments.parameters data.frame
#' @export
TLMoments.parameters.data.frame <- function(x,
                                            leftrim = attr(x, "source")$trimmings[1],
                                            rightrim = attr(x, "source")$trimmings[2],
                                            max.order = 4L,
                                            ...) {

  if (is.null(leftrim)) leftrim <- 0L
  if (is.null(rightrim)) rightrim <- 0L

  nam <- getFormulaSides(attr(x, "source")$formula)
  if (is.null(attr(x, "source")$lambdas) ||
      !identical(attr(x, "source")$trimmings, c(leftrim, rightrim)) ||
      (!is.null(attr(x, "source")$tl.order) && max(attr(x, "source")$tl.order) != max.order)) { # calculate new

    ls <- apply(x[!(names(x) %in% nam$rhs)], 1, function(xx) {
      calcTLMom(max.order, leftrim, rightrim,
                qfunc = do.call(getQ, c(x = attr(x, "distribution"), as.list(xx))))
    })
    if (max.order == 1) {
      ls <- rbind(ls, deparse.level = 0)
    }
    ls <- as.data.frame(t(ls))
    names(ls) <- paste0("L", 1:max.order)
    fac <- x[nam$rhs]

  } else { # or use old calculations
    dat <- attr(x, "source")$lambdas
    fac <- dat[nam$rhs]
    ls <- dat[!(names(dat) %in% nam$rhs)]
  }

  if (max.order > 1) {
    ratios <- t(apply(ls, 1, calcRatios))[, -1, drop = FALSE]
  } else {
    ratios <-  matrix(NA, nrow = nrow(ls), ncol = 0)
  }

  out <- list(
    lambdas = cbind(fac, as.data.frame(ls)),
    ratios = cbind(fac, as.data.frame(ratios))
  )

  do.call(returnTLMoments,c(
    list(out = out, leftrim = leftrim, rightrim = rightrim, order = 1L:max.order),
    func = "TLMoments.parameters",
    distr = attr(x, "distribution"),
    parameters = list(removeAttributes(x)),
    attr(x, "source")
  ))
}


#' @export
print.TLMoments <- function(x, ...) {
  # if ("data.frame" %in% class(x$lambdas)) {
  #   print.data.frame(x)
  #   return(invisible(x))
  # }

  tmp <- x
  attributes(tmp) <- NULL

  dim(tmp) <- dim(x)
  names(tmp) <- names(x)
  dimnames(tmp) <- dimnames(x)

  print(tmp)
  invisible(x)
}

