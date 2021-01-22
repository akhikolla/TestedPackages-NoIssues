#' @title Generalized Extreme Value distribution
#' @description Cumulative distribution function, density function, quantile function and
#' generation of random variates of the generalized extreme value distribution.
#'
#' @param x,q,p numeric vector of values, quantiles, or probabilites.
#' @param n numeric, number of random variates.
#' @param loc,scale,shape location, scale, and shape parameter of the generalized extreme value distribution. All must be of length one.
#'
#' @seealso \code{\link{pgum}}, \code{\link{pgpd}}, \code{\link{pln3}}

#' @rdname gev
#' @export
pgev <- function(q, loc = 0, scale = 1, shape = 0) {
  if (length(loc) != 1 || length(scale) != 1 || length(shape) != 1) stop("loc, scale, and shape must have length 1. ")
  if (!are.numeric(q, loc, scale, shape)) stop("q, loc, scale, and shape must be numeric. ")
  if (any(is.na(c(loc, scale, shape)))) return(rep(NA, length(q)))
  if (scale <= 0) stop("scale must be >0. ")

  y <- (q - loc)/scale
  if (abs(shape) < 1e-06) {
    exp(-exp(-y))
  } else {
    exp(-pmax(1 + shape*y, 0)^(-1/shape))
  }
}
#' @rdname gev
#' @export
dgev <- function(x, loc = 0, scale = 1, shape = 0) {
  if (length(loc) != 1 || length(scale) != 1 || length(shape) != 1) stop("loc, scale, and shape must have length 1. ")
  if (!are.numeric(x, loc, scale, shape)) stop("x, loc, scale, and shape must be numeric. ")
  if (any(is.na(c(loc, scale, shape)))) return(rep(NA, length(x)))
  if (scale <= 0) stop("scale must be >0. ")

  y <- (x - loc)/scale
  if (abs(shape) < 1e-06) {
    1/scale * exp(-y) * exp(-exp(-y))
  } else {
    1/scale * (1 + shape * y)^(-1/shape - 1) * exp(-(1 + shape*y)^(-1/shape))
  }
}
#' @rdname gev
#' @export
qgev <- function(p, loc = 0, scale = 1, shape = 0) {
  if (length(loc) != 1 || length(scale) != 1 || length(shape) != 1) stop("loc, scale, and shape must have length 1. ")
  if (!are.numeric(p, loc, scale, shape)) stop("p, loc, scale, and shape must be numeric. ")
  if (any(is.na(c(loc, scale, shape)))) return(rep(NA, length(p)))
  if (scale <= 0) stop("scale must be >0. ")

  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 1) stop("p must be between 0 and 1. ")

  if (abs(shape) < 1e-06) {
    loc - scale * log(-log(p))
  } else {
    loc + scale * ((-log(p))^(-shape) - 1)/shape
  }
}
#' @rdname gev
#' @export
rgev <- function(n, loc = 0, scale = 1, shape = 0) {
  qgev(runif(n), loc = loc, scale = scale, shape = shape)
}


#' @title Gumbel distribution
#' @description Cumulative distribution function, density function, quantile function and
#' generation of random variates of the Gumbel distribution.
#'
#' @param x,q,p numeric vector of values, quantiles, or probabilites.
#' @param n numeric, number of random variates.
#' @param loc,scale location and scale parameter of the Gumbel distribution. All must be of length one.
#'
#' @seealso \code{\link{pgev}}, \code{\link{pgpd}}, \code{\link{pln3}}

#' @rdname gum
#' @export
pgum <- function(q, loc = 0, scale = 1) {
  pgev(q, loc = loc, scale = scale, shape = 0)
}
#' @rdname gum
#' @export
dgum <- function(x, loc = 0, scale = 1) {
  dgev(x, loc = loc, scale = scale, shape = 0)
}
#' @rdname gum
#' @export
qgum <- function(p, loc = 0, scale = 1) {
  qgev(p, loc = loc, scale = scale, shape = 0)
}
#' @rdname gum
#' @export
rgum <- function(n, loc = 0, scale = 1) {
  rgev(n, loc = loc, scale = scale, shape = 0)
}


#' @title Generalized Pareto distribution
#' @description Cumulative distribution function, density function, quantile function and
#' generation of random variates of the generalized Pareto distribution.
#' @param x,q,p numeric vector of values, quantiles, or probabilites.
#' @param n numeric, number of random variates.
#' @param loc,scale,shape location, scale, and shape parameter of the generalized Pareto distribution. All must be of length one.
#'
#' @seealso \code{\link{pgev}}, \code{\link{pgum}}, \code{\link{pln3}}

#' @rdname gpd
#' @export
pgpd <- function(q, loc = 0, scale = 1, shape = 0) {
  if (length(loc) != 1 || length(scale) != 1 || length(shape) != 1) stop("loc, scale, and shape must have length 1. ")
  if (!are.numeric(q, loc, scale, shape)) stop("q, loc, scale, and shape must be numeric. ")
  if (any(is.na(c(loc, scale, shape)))) return(rep(NA, length(q)))
  if (scale <= 0) stop("scale must be >0. ")

  y <- pmax((q - loc)/scale, 0)
  if (abs(shape) < 1e-06) {
    1 - exp(-y)
  } else {
    if (shape < 0) y <- pmin(y, -1/shape)
    1 - (1 + shape*y)^(-1/shape)
  }
}
#' @rdname gpd
#' @export
dgpd <- function(x, loc = 0, scale = 1, shape = 0) {
  if (length(loc) != 1 || length(scale) != 1 || length(shape) != 1) stop("loc, scale, and shape must have length 1. ")
  if (!are.numeric(x, loc, scale, shape)) stop("x, loc, scale, and shape must be numeric. ")
  if (any(is.na(c(loc, scale, shape)))) return(rep(NA, length(x)))
  if (scale <= 0) stop("scale must be >0. ")

  y <- (x - loc)/scale
  y[y < 0] <- Inf
  if (abs(shape) < 1e-06) {
    exp(-y)/scale
  } else {
    if (shape < 0) y <- pmin(y, -1/shape)
    (1 + shape*y)^(-1/shape-1) /scale
  }
}
#' @rdname gpd
#' @export
qgpd <- function(p, loc = 0, scale = 1, shape = 0) {
  if (length(loc) != 1 || length(scale) != 1 || length(shape) != 1) stop("loc, scale, and shape must have length 1. ")
  if (!are.numeric(p, loc, scale, shape)) stop("p, loc, scale, and shape must be numeric. ")
  if (any(is.na(c(loc, scale, shape)))) return(rep(NA, length(p)))
  if (scale <= 0) stop("scale must be >0. ")

  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 1) stop("p must be between 0 and 1. ")

  if (abs(shape) < 1e-06) {
    loc - scale*log(1-p)
  } else {
    loc + scale/shape * ((1-p)^(-shape) - 1)
  }
}
#' @rdname gpd
#' @export
rgpd <- function(n, loc = 0, scale = 1, shape = 0) {
  qgpd(runif(n), loc = loc, scale = scale, shape = shape)
}


qp3 <- function(p, loc = 0, scale = 1, shape = 0) {
  if (length(loc) != 1 || length(scale) != 1 || length(shape) != 1) stop("loc, scale, and shape must have length 1. ")
  if (!are.numeric(p, loc, scale, shape)) stop("p, loc, scale, and shape must be numeric. ")
  if (any(is.na(c(loc, scale, shape)))) return(rep(NA, length(p)))
  if (scale <= 0) stop("scale must be >0. ")

  if (abs(shape) < 1e-06) { # shape = 0
    loc + scale * stats::qnorm(p)
  } else {
    if (shape > 0) {
      loc - 4 / shape^2 * abs(0.5 * scale * shape) + stats::qgamma(p, 4 / shape^2, 1 / abs(0.5 * scale * shape))
    }
    else {
      loc + 4 / shape^2 * abs(0.5 * scale * shape) - stats::qgamma(1 - p, 4 / shape^2, 1 / abs(0.5 * scale * shape))
    }
  }
}



#' @title Three-parameter lognormal distribution
#' @description Cumulative distribution function, density function, quantile function and
#' generation of random variates of the three-parameter lognormal distribution.
#' @param x,q,p numeric vector of values, quantiles, or probabilites.
#' @param n numeric, number of random variates.
#' @param loc,scale,shape location, scale, and shape parameter of the three-parameter lognormal distribution. All must be of length one.
#'
#' @seealso \code{\link{pgev}}, \code{\link{pgum}}, \code{\link{pgpd}}

#' @rdname ln3
#' @export
pln3 <- function(q, loc = 0, scale = 1, shape = 1) {
  if (length(loc) != 1 || length(scale) != 1 || length(shape) != 1) stop("loc, scale, and shape must have length 1. ")
  if (!are.numeric(q, loc, scale, shape)) stop("q, loc, scale, and shape must be numeric. ")
  if (any(is.na(c(loc, scale, shape)))) return(rep(NA, length(q)))
  #if (scale <= 0) stop("scale must be >0. ")
  if (shape <= 0) stop("shape must be >0. ")

  stats::plnorm(q - loc, meanlog = scale, sdlog = shape)
}
#' @rdname ln3
#' @export
dln3 <- function(x, loc = 0, scale = 1, shape = 1) {
  if (length(loc) != 1 || length(scale) != 1 || length(shape) != 1) stop("loc, scale, and shape must have length 1. ")
  if (!are.numeric(x, loc, scale, shape)) stop("x, loc, scale, and shape must be numeric. ")
  if (any(is.na(c(loc, scale, shape)))) return(rep(NA, length(x)))
  #if (scale <= 0) stop("scale must be >0. ")
  if (shape <= 0) stop("shape must be >0. ")

  stats::dlnorm(x - loc, meanlog = scale, sdlog = shape)
}
#' @rdname ln3
#' @export
qln3 <- function(p, loc = 0, scale = 1, shape = 1) {
  if (length(loc) != 1 || length(scale) != 1 || length(shape) != 1) stop("loc, scale, and shape must have length 1. ")
  if (!are.numeric(p, loc, scale, shape)) stop("p, loc, scale, and shape must be numeric. ")
  if (any(is.na(c(loc, scale, shape)))) return(rep(NA, length(p)))
  #if (scale <= 0) stop("scale must be >0. ")
  if (shape <= 0) stop("shape must be >0. ")

  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 1) stop("p must be between 0 and 1. ")

  stats::qlnorm(p, meanlog = scale, sdlog = shape) + loc
}
#' @rdname ln3
#' @export
rln3 <- function(n, loc = 0, scale = 1, shape = 1) {
  if (length(loc) != 1 || length(scale) != 1 || length(shape) != 1) stop("loc, scale, and shape must have length 1. ")
  if (!are.numeric(loc, scale, shape)) stop("loc, scale, and shape must be numeric. ")
  if (!are.integer.like(n)) stop("n must be integer-like. ")
  #if (scale <= 0) stop("scale must be >0. ")
  if (shape <= 0) stop("shape must be >0. ")

  stats::rlnorm(n, meanlog = scale, sdlog = shape) + loc
}


