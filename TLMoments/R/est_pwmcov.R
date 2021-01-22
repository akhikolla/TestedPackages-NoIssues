#' @title
#' Estimate the covariance matrix of PWM estimations
#' @description
#' Internal function. Use \link{est_cov}. Description not done yet.
#'
#' @param x numeric vector or matrix of data.
#' @param order numeric vector giving the orders that are returned.
#' @param distr character of length 1 which indicates a distribution if a
#' parametric assumption should be used.
#' @param distr.trim integer vector of length 2 indicating the trimming used
#' to calculate parameters if a parametric assumption is used (i.e. \code{distr} is set).
#'
#' @return numeric matrix
#'
#' @examples
#' ### Numeric vectors
#' x <- rgev(500, shape = .2)
#' est_pwmcov(x)
#' est_pwmcov(x, distr = "gev")
#'
#' ### Numeric matrices
#' x <- matrix(rgev(600, shape = .2), nc = 3)
#' est_pwmcov(x, order = 0:2)
#' est_pwmcov(x, order = 0:2, distr = "gev")
#'
#' @rdname est_pwmcov
#' @export
est_pwmcov <- function(x, order = 0:3, distr = NULL, distr.trim = c(0, 0)) {
  if (inherits(x, c("PWMs", "TLMoments", "parameters", "quantiles")))
    stop("est_pwmcov is only for data vectors or matrices. ")

  if (is.double(order)) order <- as.integer(order)
  if (!is.integer(order)) stop("order has to be integer type!")

  if (!is.null(distr) && distr != "gev")
    stop("distr has to be NULL (nonparametric) or \"gev\"")
  if (!is.null(distr) && (!are.integer.like(distr.trim[1], distr.trim[2]) | length(distr.trim) != 2))
    stop("If distr is not NULL, distr.trim has to be a integer-like vector of length 2. ")

  UseMethod("est_pwmcov")
}

#' @rdname est_pwmcov
#' @method est_pwmcov numeric
#' @export
est_pwmcov.numeric <- function(x, order = 0:3, distr = NULL, distr.trim = c(0, 0)) {

  if (is.null(distr)) {
    # nonparametric estimation
    r <- stats::cov(pseudo(x, order))
  } else {
    # parametric estimation
    # scale and shape are estimated using L-moments.
    p <- parameters(TLMoments(x, leftrim = distr.trim[1], rightrim = distr.trim[2], na.rm = TRUE), distr)

    if (p["shape"] >= .5) {
      warning("estimated shape parameter is too high. ")
      r <- matrix(Inf, nrow = length(order), ncol = length(order))
    } else {
      r <- parametricPWMCov(distr, order, scale = p["scale"], shape = p["shape"])
    }
  }

  rownames(r) <- colnames(r) <- buildNames("beta", order)
  out <- r / sum(!is.na(x))
  attr(out, "distribution") <- distr
  out
}

#' @rdname est_pwmcov
#' @method est_pwmcov matrix
#' @export
est_pwmcov.matrix <- function(x, order = 0:3, distr = NULL, distr.trim = c(0, 0)) {

  J <- ncol(x) # Anzahl Pegel
  # Berechnung von Pseudobeobachtungen fuer jede Dimension
  pseudos <- lapply(1L:J, function(i) pseudo(x[, i], order))
  # Berechnung der Sigma-Matrix (ohne Korrekturvorfaktor)
  sigma <- stats::cov(do.call("cbind", pseudos), use = "pairwise.complete.obs")
  rownames(sigma) <- colnames(sigma) <- buildNames("beta", order, 1:J)

  if (!is.null(distr)) {
    # parametric estimation, overwrites non-diagonal block matrices
    # scale and shape are estimated using L-moments.
    xl <- lapply(1:ncol(x), function(i) x[, i])
    p <- parameters(TLMoments(xl, leftrim = distr.trim[1], rightrim = distr.trim[2], na.rm = TRUE), distr)
    r <- lapply(p, function(p) {
      parametricPWMCov(distr, order, scale = p["scale"], shape = p["shape"])
    })
    sigma <- blockdiag_list(r, sigma)
  }

  # Calculation correction factor
  rs <- apply(x, 2, function(y) sum(!is.na(y))/length(y))
  vorf <- do.call(rbind, lapply(1:J, function(i) {
    do.call(cbind, lapply(1:J, function(j) {
      matrix(min(c(rs[i], rs[j]))/(rs[i] * rs[j]), nrow = length(order), ncol = length(order))
    }))
  }))

  # Ausgeben
  out <- vorf*sigma / max(apply(x, 2, function(y) sum(!is.na(y))))
  attr(out, "distribution") <- distr
  out
}

