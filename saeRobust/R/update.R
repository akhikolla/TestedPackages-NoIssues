#' Update a fitted object
#'
#' This is a method which can be used to update a \link{rfh} result object and
#' refit it. The fitted parameter values from the current object are used as
#' starting values, then \link{update.default} is called.
#'
#' @param object (rfh) an object fitted by \link{rfh}
#' @param ... arguments passed to \link{update.default}
#' @param formula see \link{update.formula}
#' @param where (environment) should not be specified by the user
#'
#' @include NAMESPACE.R
#'
#' @rdname update
#' @export
update(object ~ rfh, formula, ..., where = parent.frame(2)) %m% {

  object$call[c("x0Coef", "x0Var", "x0Re")] <-
    object[c("coefficients", "variance", "re")]

  fieldNames <- c("k", "tol", "maxIter", "maxIterParam", "maxIterRe")
  object$call[fieldNames] <- object[fieldNames]

  # continuing with update.default because of nse problems
  call <- getCall(object)
  extras <- list(...)

  if (!missing(formula) && inherits(formula, "formula"))
    call$formula <- update(formula(object), formula)

  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }

  eval(call, where)

}

#' @rdname update
#' @export
update(object ~ fitrfh, ...) %m% {
  # the first class should be the fitting function. If it is a rfh object this
  # method should never be called
  fun <- class(object)[1] # yeah...

  object[c("x0Coef", "x0Var", "x0Re")] <-
    object[c("coefficients", "variance", "re")]

  args <- list(...)
  object[names(args)] <- args

  do.call(get(fun, mode = "function"), object)

}
