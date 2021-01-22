#' @title
#' Estimate the covariance matrix of quantile estimations
#' @description
#' Internal function. Use \link{est_cov}. Description not done yet.
#'
#' @param x numeric vector or matrix containing data.
#' @param distr character of length 1 giving the distribution if parametric assumption should be used.
#' @param p quantile levels from which the covariance should be calculated.
#' @param leftrim,rightrim lower and upper trimming parameter used for parameter calculation, have to be non-negative integers.
#' @param np.cov boolean, if TRUE no parametric assumptions are used to calculate the covariance matrix (default FALSE).
#' @param reg.weights numeric vector of weights for regionalized TLMoments.
#' @param set.n hypothetical data length n if theoretical values are given.
#' @param ... additional arguments.
#'
#' @return numeric matrix
#'
#' @examples
#' ### Numeric vectors
#' x <- rgev(500, shape = .2)
#'
#' quantiles(parameters(TLMoments(x), "gev"), c(.9, .95, .99))
#' est_quancov(x, "gev", c(.9, .95, .99), 0, 0)
#' #cov(t(replicate(5000,
#' #  quantiles(parameters(TLMoments(rgev(500, shape = .2)), "gev"), c(.9, .95, .99))
#' #)))
#'
#' quantiles(parameters(TLMoments(x, rightrim = 1), "gev"), c(.9, .95, .99))
#' est_quancov(x, "gev", c(.9, .95, .99), 0, 1)
#' #cov(t(replicate(5000,
#' #  quantiles(
#' #    parameters(TLMoments(rgev(500, shape = .2), rightrim = 1), "gev"),
#' #    c(.9, .95, .99)
#' #  )
#' #)))
#'
#' ### Numeric matrices
#' x <- matrix(rgev(600, shape = .2), nc = 3)
#'
#' quantiles(parameters(TLMoments(x), "gev"), c(.9, .95, .99))
#' est_quancov(x, "gev", c(.9, .95, .99), 0, 0)
#'
#' est_quancov(x, "gev", .9, 0, 0)
#' #cov(t(replicate(5000,
#' # quantiles(
#' #   parameters(TLMoments(matrix(rgev(600, shape = .2), nc = 3)),
#' #  "gev"), .9)
#' #  )
#' #))
#'
#' ### quantiles object
#' q <- quantiles(as.parameters(loc = 3, scale = 2, shape = .4, distr = "gev"), c(.9, .99))
#' est_quancov(q)
#' est_quancov(q, leftrim = 0, rightrim = 0)
#' est_quancov(q, leftrim = 0, rightrim = 0, set.n = 10)
#'
#' @rdname est_quancov
#' @export
est_quancov <- function(x,
                        distr = "",
                        p = NULL,
                        leftrim = 0L,
                        rightrim = 0L,
                        ...) {
  if (inherits(x, "quantiles")) {
    distr <- attr(x, "distribution")
    p <- attr(x, "p")
    if (!any(is.null(attr(x, "source")$trimming))) {
      leftrim <- attr(x, "source")$trimming[1]
      rightrim <- attr(x, "source")$trimming[2]
    }
  }

  # TODO: General error catches
  if (is.na(leftrim) | is.na(rightrim))
    stop("leftrim and/or right must be given. ")

  if (!are.integer.like(leftrim, rightrim))
    stop("leftrim and rightrim have to be integer types!")

  if (!(distr %in% c("gev", "gum", "gpd", "ln3")))
    stop("distr has to be \"gev\", \"gum\", \"gpd\", or \"ln3\". ")

  if (!(distr %in% c("gev")))
    stop("only GEV for now")

  if (!is.numeric(p))
    stop("Argument p must be numeric vector. ")

  UseMethod("est_quancov")
}

#' @rdname est_quancov
#' @method est_quancov numeric
#' @export
est_quancov.numeric <- function(x,
                                distr,
                                p,
                                leftrim = 0L,
                                rightrim = 0L,
                                np.cov = FALSE,
                                ...) {

  if (np.cov) {
    paramcov <- est_paramcov(x, distr, leftrim = leftrim, rightrim = rightrim, np.cov = TRUE)
  } else {
    paramcov <- est_paramcov(x, distr, leftrim = leftrim, rightrim = rightrim)
  }
  param <- parameters(TLMoments(x, leftrim = leftrim, rightrim = rightrim, na.rm = TRUE), distr)
  A <- CovParamtoQuan(distr, param, p)
  r <- A %*% paramcov %*% t(A)
  rownames(r) <- colnames(r) <- buildNames("q", p)

  attr(r, "distribution") <- distr
  attr(r, "trimmings") <- c(leftrim, rightrim)
  attr(r, "n") <- length(x)
  r
}

#' @rdname est_quancov
#' @method est_quancov matrix
#' @export
est_quancov.matrix <- function(x,
                               distr,
                               p,
                               leftrim = 0L,
                               rightrim = 0L,
                               np.cov = FALSE,
                               reg.weights = NULL,
                               ...) {

  if (np.cov) {
    paramcov <- est_paramcov(x, distr, leftrim, rightrim, np.cov = TRUE, reg.weights = reg.weights)
  } else {
    paramcov <- est_paramcov(x, distr, leftrim, rightrim, reg.weights = reg.weights)
  }

  if (is.null(reg.weights)) {
    param <- parameters(TLMoments(x, leftrim = leftrim, rightrim = rightrim, na.rm = TRUE), distr)
    A <- lapply(1:ncol(param), function(i) CovParamtoQuan(distr, param[, i], p))
    A <- blockdiag_list(A)
    r <- A %*% paramcov %*% t(A)
    rownames(r) <- colnames(r) <- buildNames("q", p, 1:ncol(x))
  } else {
    param <- parameters(regionalize(TLMoments(x, leftrim = leftrim, rightrim = rightrim, na.rm = TRUE), reg.weights), distr)
    A <- CovParamtoQuan(distr, param, p)
    r <- A %*% paramcov %*% t(A)
    rownames(r) <- colnames(r) <- buildNames("q", p)
  }

  attr(r, "distribution") <- distr
  attr(r, "trimmings") <- c(leftrim, rightrim)
  attr(r, "n") <- apply(x, 2, function(y) sum(!is.na(y)))
  r
}

#' @rdname est_quancov
#' @method est_quancov quantiles
#' @export
est_quancov.quantiles <- function(x,
                                  distr = attr(x, "distribution"),
                                  p = attr(x, "p"),
                                  leftrim = attr(x, "source")$trimmings[1],
                                  rightrim = attr(x, "source")$trimmings[2],
                                  set.n = NA,
                                  ...) {
  if (!inherits(x, "numeric"))
    stop("x must be a numeric quantiles-object. ")
  if (!any(attr(x, "source")$func %in% c("as.PWMs", "as.TLMoments", "as.parameters")))
    stop("est_quancov.quantiles only for theoretical values. ")

  if (are.integer.like(leftrim) && is.null(rightrim)) {
    rightrim <- 0
  }
  if (are.integer.like(rightrim) && is.null(leftrim)) {
    leftrim <- 0
  }
  if (any(is.null(c(leftrim, rightrim)))) {
    #warning("No trimmings found, set to leftrim = 0, rightrim = 0. ")
    leftrim <- rightrim <- 0
  }

  if (is.na(set.n) | !is.numeric(set.n)) {
    #warning("Missing or invalid set.n argument. Giving results for n = 100. ")
    n <- 100
  } else { n <- set.n }

  param <- as.parameters(attr(x, "source")$parameters, distr = distr)
  paramcov <- est_paramcov(param, leftrim = leftrim, rightrim = rightrim, set.n = n)
  A <- CovParamtoQuan(distr, param, p)
  r <- A %*% paramcov %*% t(A)
  rownames(r) <- colnames(r) <- buildNames("q", p)

  attr(r, "distribution") <- distr
  attr(r, "trimmings") <- c(leftrim, rightrim)
  attr(r, "n") <- n
  r
}
