#' @title
#' Covariance matrix of PWMs, TLMoments, parameters, or quantiles
#' @description
#' Calculation of the empirical or theoretical covariance matrix of objects
#' of the classes \code{PWMs}, \code{TLMoments}, \code{parameters}, or \code{quantiles}.
#'
#' @param x object of \code{PWMs}, \code{TLMoments}, \code{parameters},
#' or \code{quantiles} constructed using the same-named functions.
#' @param select numeric oder character vector specifying a subset of the covariance matrix. If
#' not specified the full covariance matrix is returned.
#' @param ... additional arguments given to the sub-functions: \code{distr} and \code{np.cov} (see details).
#'
#' @details
#' Covariance matrices of \code{PWMs} and \code{TLMoments} are calculated without parametric
#' assumption by default. Covariance matrices of \code{parameters} and \code{quantiles} use
#' parametric assumption based on their stored distribution attribute (only GEV at the moment).
#' Parametric (GEV) calculation can be enforced by specifying \code{distr=\"gev\"}, non-parametric
#' calculation by using \code{np.cov=TRUE}.
#'
#' @return numeric matrix (if \code{x} is of class \code{PWMs}, \code{parameters}, or
#'          \code{quantiles}) or a list of two matrices (\code{lambdas} and \code{ratios}, if
#'          \code{x} is of class \code{TLMoments}).
#'
#' @seealso \code{\link{PWMs}}, \code{\link{TLMoments}}, \code{\link{parameters}}, \code{\link{quantiles}}
#' @examples
#' ### 1: PWMs:
#'
#' xvec <- rgev(100, shape = .1)
#' xmat <- cbind(rgev(100, shape = .1), rgev(100, shape = .3))
#'
#' # Covariance estimation of PWMs normally without parametric assumption:
#' est_cov(PWMs(xvec))
#' est_cov(PWMs(xvec), select = 0:1)
#' est_cov(PWMs(xmat))
#' est_cov(PWMs(xmat), select = 3)
#' est_cov(PWMs(xmat[, 1, drop = FALSE]), select = 2:3)
#'
#' # Parametric assumptions (only GEV by now) can be used:
#' est_cov(PWMs(xvec), distr = "gev")
#' est_cov(PWMs(xvec), distr = "gev", select = c(1, 3))
#'
#' \dontrun{
#' cov(t(replicate(100000,
#'   as.vector(PWMs(cbind(rgev(100, shape = .1), rgev(100, shape = .3)), max.order = 1)))
#' ))
#' }
#'
#'
#' ### 2. TLMoments:
#'
#' xvec <- rgev(100, shape = .1)
#' xmat <- cbind(rgev(100, shape = .1), rgev(100, shape = .3))
#'
#' # Covariance estimation of TLMoments normally without parametric assumption:
#' est_cov(TLMoments(xvec))
#' est_cov(TLMoments(xvec, rightrim = 1))
#' est_cov(TLMoments(xvec), select = 3:4)
#'
#' # Parametric assumptions (only GEV by now) can be used:
#' est_cov(TLMoments(xvec), distr = "gev")
#'
#' # Matrix inputs
#' est_cov(TLMoments(xmat))
#' est_cov(TLMoments(xmat), select = 3:4)
#' est_cov(TLMoments(xmat[, 1, drop = FALSE]), select = 3:4)
#'
#' # Covariance of theoretical TLMoments only with parametric assumption:
#' est_cov(as.TLMoments(c(14, 4, 1)), distr = "gev", set.n = 100)
#' est_cov(as.TLMoments(c(14, 4, 1), rightrim = 1), distr = "gev", set.n = 100)
#'
#' # Regionalized TLMoments
#' est_cov(regionalize(TLMoments(xmat), c(.75, .25)))
#' est_cov(regionalize(TLMoments(xmat), c(.75, .25)), distr = "gev", select = 3:4)
#'
#'
#' ### 3. Parameters:
#'
#' xvec <- rgev(100, shape = .1)
#' xmat <- cbind(rgev(100, shape = .1), rgev(100, shape = .3))
#'
#' # Covariance estimation of parameters normally with parametric assumption:
#' est_cov(parameters(TLMoments(xvec), "gev"))
#' est_cov(parameters(TLMoments(xvec, rightrim = 1), "gev"))
#' est_cov(parameters(TLMoments(xvec, rightrim = 1), "gev"), select = c("scale", "shape"))
#'
#' # A nonparametric estimation can be enforced with np.cov:
#' est_cov(parameters(TLMoments(xvec), "gev"), np.cov = TRUE)
#' est_cov(parameters(TLMoments(xvec, rightrim = 1), "gev"), np.cov = TRUE)
#'
#' # Matrix inputs
#' est_cov(parameters(TLMoments(xmat), "gev"))
#' est_cov(parameters(TLMoments(xmat), "gev"), select = "shape")
#' est_cov(parameters(TLMoments(xmat[, 1]), "gev"), select = "shape")
#'
#' # Theoretical values (leftrim and/or rightrim have to be specified)
#' para <- as.parameters(loc = 10, scale = 5, shape = .2, distr = "gev")
#' est_cov(para, set.n = 100)
#' est_cov(para, rightrim = 1, set.n = 100)
#'
#' \dontrun{
#' var(t(replicate(10000, parameters(TLMoments(rgev(100, 10, 5, .2)), "gev"))))
#' }
#' \dontrun{
#' var(t(replicate(10000, parameters(TLMoments(rgev(100, 10, 5, .2), rightrim = 1), "gev"))))
#' }
#'
#' # Parameter estimates from regionalized TLMoments:
#' est_cov(parameters(regionalize(TLMoments(xmat), c(.75, .25)), "gev"))
#'
#'
#' ### 4. Quantiles:
#'
#' xvec <- rgev(100, shape = .2)
#' xmat <- cbind(rgev(100, shape = .1), rgev(100, shape = .3))
#'
#' # Covariance estimation of parameters normally with parametric assumption:
#' q <- quantiles(parameters(TLMoments(xvec), "gev"), c(.9, .95, .99))
#' est_cov(q)
#' est_cov(q, select = c("0.9", "0.99"))
#' est_cov(q, select = .95)
#'
#' # A nonparametric estimation can be enforced with np.cov:
#' est_cov(q, np.cov = TRUE)
#'
#' # Matrix inputs
#' param <- parameters(TLMoments(xmat, 0, 1), "gev")
#' q <- quantiles(param, c(.9, .95, .99))
#' est_cov(q)
#' est_cov(q, select = .99)
#' param <- parameters(TLMoments(xmat[, 1, drop = FALSE], 0, 1), "gev")
#' q <- quantiles(param, c(.9, .95, .99))
#' est_cov(q, select = .99)
#'
#' # Theoretical values
#' q <- quantiles(as.parameters(loc = 10, scale = 5, shape = .3, distr = "gev"), c(.9, .99))
#' est_cov(q)
#' est_cov(q, leftrim = 0, rightrim = 1)
#' est_cov(q, leftrim = 0, rightrim = 1, set.n = 100)
#'
#' # Quantile estimates from regionalized TLMoments:
#' param <- parameters(regionalize(TLMoments(xmat), c(.75, .25)), "gev")
#' est_cov(quantiles(param, c(.9, .99)))
#'
#'
#' @rdname est_cov
#' @export
est_cov <- function(x, ...) {
  # TODO: General error catches
  if (!inherits(x, c("PWMs", "TLMoments", "parameters", "quantiles")))
    stop("x must be of class PWMs, TLMoments, parameters, or quantiles. ")

  UseMethod("est_cov")
}

#' @rdname est_cov
#' @method est_cov PWMs
#' @export
est_cov.PWMs <- function(x,
                         select = attr(x, "order"),
                         ...) {

  # TODO: General error catches
  if (!inherits(x, c("numeric", "matrix")))
    stop("To date, only numeric and matrix types of PWMs are supported. ")


  ret <- est_pwmcov(attr(x, "source")$data, order = attr(x, "order"), ...)
  i <- grep(paste0("beta(", paste(select, collapse = "|"), "){1}(_[0-9]*)?$"), rownames(ret))
  out <- ret[i, i, drop = FALSE]

  attr(out, "distribution") <- attr(ret, "distribution")
  class(out) <- "est_cov"
  out
}

#' @rdname est_cov
#' @method est_cov TLMoments
#' @export
est_cov.TLMoments <- function(x,
                              select = attr(x, "order"),
                              ...) {

  # TODO: General error catches
  if (!inherits(x$lambdas, c("numeric", "matrix")))
    stop("To date, only numeric and matrix types of TLMoments are supported. ")

  if (any(attr(x, "source")$func %in% c("as.PWMs", "as.TLMoments", "as.parameters"))) { # theoretical values
    ret <- est_tlmcov(x, ...)
  } else { # empirical values
    ret <- est_tlmcov(attr(x, "source")$data,
               leftrim = attr(x, "leftrim"),
               rightrim = attr(x, "rightrim"),
               order = attr(x, "order"),
               reg.weights = attr(x, "source")$reg.weights,
               ...)
  }

  # est_tlmcov can return list or single matrix (by setting lambda.cov or ratio.cov to FALSE)
  if (is.list(ret)) {
    i1 <- grep(paste0("^L(", paste(select, collapse = "|"), "){1}(_[0-9]*)?$"), rownames(ret$lambdas))
    i2 <- grep(paste0("^T(", paste(select, collapse = "|"), "){1}(_[0-9]*)?$"), rownames(ret$ratios))
    out <- list(lambdas = ret$lambdas[i1, i1, drop = FALSE],
                ratios = ret$ratios[i2, i2, drop = FALSE])
  } else {
    i <- grep(paste0("^(L|T)(", paste(select, collapse = "|"), "){1}(_[0-9]*)?$"), rownames(ret))
    out <- ret[i, i, drop = FALSE]
  }
  attr(out, "distribution") <- attr(ret, "distribution")
  attr(out, "trimmings") <- attr(ret, "trimmings")
  attr(out, "n") <- attr(ret, "n")
  class(out) <- "est_cov"
  out
}

#' @rdname est_cov
#' @method est_cov parameters
#' @export
est_cov.parameters <- function(x,
                               select = c("loc", "scale", "shape"),
                               ...) {

  # TODO: General error catches
  if (!inherits(x, c("numeric", "matrix")))
    stop("To date, only numeric and matrix types of parameters are supported. ")

  if (any(attr(x, "source")$func %in% c("as.PWMs", "as.TLMoments", "as.parameters"))) {  # theoretical values
    # if (!are.integer.like(attr(x, "source")$trimmings)) {
    #   attr(x, "source")$trimmings <- c(leftrim, rightrim)
    #   #warning("Calculation based on TL(0,0)-moments. Overwrite with \"leftrim\" and \"rightrim\". ")
    # }
    ret <- est_paramcov(x, ...)
  } else { # empirical values
    ret <- est_paramcov(attr(x, "source")$data,
                 distr = attr(x, "distribution"),
                 leftrim = attr(x, "source")$trimmings[1],
                 rightrim = attr(x, "source")$trimmings[2],
                 reg.weights = attr(x, "source")$reg.weights,
                 ...)
  }

  i <- grep(paste0(paste(select, collapse = "|"), "(_[0-9]*)?$"), rownames(ret))
  out <- ret[i, i, drop = FALSE]
  attr(out, "distribution") <- attr(ret, "distribution")
  attr(out, "trimmings") <- attr(ret, "trimmings")
  attr(out, "n") <- attr(ret, "n")
  class(out) <- "est_cov"
  out
}

#' @rdname est_cov
#' @method est_cov quantiles
#' @export
est_cov.quantiles <- function(x,
                              select = attr(x, "p"),
                              ...) {

  # TODO: General error catches
  if (!inherits(x, c("numeric", "matrix")))
    stop("To date, only numeric and matrix types of quantiles are supported. ")

  if (any(attr(x, "source")$func %in% c("as.PWMs", "as.TLMoments", "as.parameters"))) {  # theoretical values
    # if (!are.integer.like(attr(x, "source")$trimmings)) {
    #   attr(x, "source")$trimmings <- c(leftrim, rightrim)
    #   #warning("Calculation based on TL(0,0)-moments. Overwrite with \"leftrim\" and \"rightrim\". ")
    # }
    ret <- est_quancov(x, ...)
  } else { # empirical values
    ret <- est_quancov(attr(x, "source")$data,
                distr = attr(x, "distribution"),
                p = attr(x, "p"),
                leftrim = attr(x, "source")$trimmings[1],
                rightrim = attr(x, "source")$trimmings[2],
                reg.weights = attr(x, "source")$reg.weights,
                ...)
  }

  dims <- sub("^q([0-9]\\.[0-9]*)(_[0-9]*)?$", "\\1", x = rownames(ret))
  i <- dims %in% as.character(select)
  out <- ret[i, i, drop = FALSE]

  attr(out, "distribution") <- attr(ret, "distribution")
  attr(out, "trimmings") <- attr(ret, "trimmings")
  attr(out, "n") <- attr(ret, "n")
  class(out) <- "est_cov"
  out
}

#' @export
print.est_cov <- function(x, ...) {
  tmp <- x
  attributes(tmp) <- NULL
  dim(tmp) <- dim(x)
  names(tmp) <- names(x)
  dimnames(tmp) <- dimnames(x)

  print(tmp)
  invisible(x)
}
