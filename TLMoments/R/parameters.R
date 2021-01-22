#' @title
#' Converting TL-moments to distribution parameters
#' @description
#' Converts TL-moments (or PWMs) to distribution parameters. By now, conversions for gev, gumbel, gpd, and ln3
#' are available. Important trimming options are calculated through known formulas (see references for
#' some of them), other options are calculated through a numerical optimization.
#'
#' @param x object returned by TLMoments (or PWMs, in which case TLMoments with no trimming is used).
#' @param distr character object defining the distribution. Supported types are
#' "gev", "gum", "gpd", and "ln3".
#' @param ... additional arguments.
#'
#' @return numeric vector, matrix, list, or data.frame of parameter estimates with
#' class \code{parameters}. The object contains the following attributes: \itemize{
#'  \item \code{distribution}: a character indicating the used distribution
#'  \item \code{source}: a list with background information (used function, data, n, formula, trimmings; mainly for internal purposes)
#' }
#' The attributes are hidden in the print-function for a clearer presentation.
#'
#' @references Elamir, E. A. H. (2010). Optimal choices for trimming in trimmed L-moment method. Applied Mathematical Sciences, 4(58), 2881-2890.
#' @references Fischer, S., Fried, R., & Schumann, A. (2015). Examination for robustness of parametric estimators for flood statistics in the context of extraordinary extreme events. Hydrology and Earth System Sciences Discussions, 12, 8553-8576.
#' @references Hosking, J. R. (1990). L-moments: analysis and estimation of distributions using linear combinations of order statistics. Journal of the Royal Statistical Society. Series B (Methodological), 105-124.
#' @references Hosking, J. R. M. (2007). Some theory and practical uses of trimmed L-moments. Journal of Statistical Planning and Inference, 137(9), 3024-3039.
#' @seealso \code{\link{PWMs}}, \code{\link{TLMoments}}, \code{\link{quantiles}}, \code{\link{summary.parameters}}, \code{\link{as.parameters}}. Built-in distributions: \code{\link{pgev}}, \code{\link{pgum}}, \code{\link{pgpd}}, \code{\link{pln3}}.
#'
#' @examples
#' xmat <- matrix(rgev(100, shape = .2), nc = 4)
#' xvec <- xmat[, 3]
#' xlist <- lapply(1L:ncol(xmat), function(i) xmat[, i])
#' xdat <- data.frame(
#'  station = rep(letters[1:2], each = 50),
#'  season = rep(c("S", "W"), 50),
#'  hq = as.vector(xmat)
#' )
#'
#' # TLMoments-objects or PWMs-objects can be used. However, in case of PWMs
#' # simply the TLMoments(., leftrim = 0, rightrim = 0)-variant is used.
#'
#' parameters(PWMs(xvec), "gev")
#' tlm <- TLMoments(xvec, leftrim = 0, rightrim = 0)
#' parameters(tlm, "gev")
#'
#' tlm <- TLMoments(xmat, leftrim = 1, rightrim = 1)
#' parameters(tlm, "gum")
#'
#' tlm <- TLMoments(xlist)
#' parameters(tlm, "gpd")
#'
#' tlm <- TLMoments(xdat, hq ~ station, leftrim = 0, rightrim = 2)
#' parameters(tlm, "gev")
#'
#' tlm <- TLMoments(xdat, hq ~ station + season, leftrim = 0, rightrim = 2)
#' parameters(tlm, "gev")
#'
#'
#' # If no explicit formula is implemented, it is tried to calculate
#' # parameters numerically. The attribute source$param.computation.method
#' # indicates if this is the case.
#'
#' tlm <- TLMoments(rgum(200, loc = 5, scale = 2), leftrim = 1, rightrim = 4)
#' parameters(tlm, "gum")
#'
#' tlm <- TLMoments(rgev(200, loc = 10, scale = 5, shape = .4), leftrim = 2, rightrim = 2)
#' parameters(tlm, "gev")
#'
#' tlm <- TLMoments(rln3(200, loc = 3, scale = 1.5, shape = 2), leftrim = 0, rightrim = 1)
#' parameters(tlm, "ln3")
#'
#' # Numerical calculation is A LOT slower:
#' \dontrun{
#' system.time(replicate(500,
#'   parameters(TLMoments(rgum(100, loc = 5, scale = 2), 1, 1), "gum")
#' ))[3]
#' system.time(replicate(500,
#'   parameters(TLMoments(rgum(100, loc = 5, scale = 2), 1, 2), "gum")
#' ))[3]
#' }
#'
#' # Using magrittr
#' library(magrittr)
#'
#' TLMoments(rgpd(500, loc = 10, scale = 3, shape = .3), rightrim = 0) %>%
#'  parameters("gpd")
#'
#' TLMoments(rgpd(500, loc = 10, scale = 3, shape = .3), rightrim = 0) %>%
#'  parameters("gpd", u = 10)
#'
#' TLMoments(rgpd(500, loc = 10, scale = 3, shape = .3), rightrim = 1) %>%
#'  parameters("gpd")
#'
#' TLMoments(rgpd(500, loc = 10, scale = 3, shape = .3), rightrim = 2) %>%
#'  parameters("gpd")
#'
#' @export
parameters <- function(x, distr, ...) {
  if (!inherits(x, c("TLMoments", "PWMs"))) stop("tlm has to be of class TLMoments or PWMs")

  valid.distr <- c("gev", "gum", "gpd", "ln3")
  if (!(distr %in% valid.distr))
    stop(paste0("distr has to be one of ", paste(valid.distr, collapse = ", "), "."))

  # Check for sufficient data:
  if (distr == "gev") {
    if (!all(c(1, 2, 3) %in% attr(x, "order"))) stop("Parameter calculation of \"gev\" needs at least the first three TL-moments")
  } else if (distr == "gum") {
    if (!all(c(1, 2) %in% attr(x, "order"))) stop("Parameter calculation of \"gum\" needs at least the first two TL-moments")
  } else if (distr == "gpd" && "u" %in% names(list(...))) {
    if (!all(c(1, 2) %in% attr(x, "order"))) stop("Parameter calculation of \"gpd\" needs at least the first two TL-moments")
  } else if (distr == "gpd") {
    if (!all(c(1, 2, 3) %in% attr(x, "order"))) stop("Parameter calculation of \"gpd\" needs at least the first three TL-moments or a specified threshold value (u = ...)")
  } else if (distr == "ln3") {
    if (!all(c(1, 2, 3) %in% attr(x, "order"))) stop("Parameter calculation of \"ln3\" needs at least the first three TL-moments")
  }

  UseMethod("parameters")
}

#' @title returnParameters
#' @description Sets attributions to parameters objects and returns them. This function is for internal use.
#' @param out -
#' @param distribution -
#' @param ... -
#' @return An object of class parameters.
returnParameters <- function(out, distribution, ...) {

  class <- class(out)
  args <- list(...)

  # If no func attribute is set, set to
  if (!exists("func", args)) args$func <- "parameters"

  # If more than one func attributes are given, concatenate them
  if (sum(names(args) == "func") >= 2) {
    newfunc <- as.vector(unlist(args[names(args) == "func"]))
    args$func <- NULL
    args$func <- newfunc
  }

  # Attributes of parameters
  # distribution
  # source: func
  #         data (if calculated)
  #         input (if not calculated)
  #         n (if calculated)
  #         formula (if data is data.frame)
  #         trimmings (if coming from TLMoments)
  #         lambdas (if coming from TLMoments)
  #         tl.order (if coming from TLMoments)
  # class: "parameters", class

  attr(out, "distribution") <- distribution
  attr(out, "source") <- args
  class(out) <- c("parameters", class)

  out
}

#' @describeIn parameters parameters for PWMs-object
#' @method parameters PWMs
#' @export
parameters.PWMs <- function(x, distr, ...) {
  parameters.TLMoments(TLMoments(x), distr = distr)
}

#' @describeIn parameters parameters for TLMoments-object
#' @method parameters TLMoments
#' @export
parameters.TLMoments <- function(x, distr, ...) {
  UseMethod("parameters.TLMoments", x$lambdas)
}

#' @method parameters.TLMoments numeric
#' @export
parameters.TLMoments.numeric <- function(x, distr, ...) {
  leftrim <- attr(x, "leftrim")
  rightrim <- attr(x, "rightrim")

  out <- param.est(x$lambdas, leftrim, rightrim, distr, ...)
  param.computation.method <- attr(out, "param.computation.method")
  attr(out, "param.computation.method") <- NULL

  do.call(returnParameters, c(
    list(out = out, distribution = distr),
    func = "parameters",
    lambdas = list(x$lambdas),
    trimmings = list(c(attr(x, "leftrim"), attr(x, "rightrim"))),
    tl.order = list(attr(x, "order")),
    param.computation.method = param.computation.method,
    attr(x, "source")
  ))
}

#' @method parameters.TLMoments matrix
#' @export
parameters.TLMoments.matrix <- function(x, distr, ...) {
  leftrim <- attr(x, "leftrim")
  rightrim <- attr(x, "rightrim")

  out <- lapply(1L:ncol(x$lambdas), function(i) {
    param.est(x$lambdas[, i], leftrim, rightrim, distr, ...)
  })
  param.computation.method <- attr(out[[1]], "param.computation.method")
  out <- simplify2array(out)

  do.call(returnParameters, c(
    list(out = out, distribution = distr),
    func = "parameters",
    lambdas = list(x$lambdas),
    trimmings = list(c(attr(x, "leftrim"), attr(x, "rightrim"))),
    tl.order = list(attr(x, "order")),
    param.computation.method = param.computation.method,
    attr(x, "source")
  ))
}

#' @method parameters.TLMoments list
#' @export
parameters.TLMoments.list <- function(x, distr, ...) {
  leftrim <- attr(x, "leftrim")
  rightrim <- attr(x, "rightrim")

  out <- lapply(seq_along(x$lambdas), function(i) {
    param.est(x$lambdas[[i]], leftrim, rightrim, distr, ...)
  })
  param.computation.method <- attr(out[[1]], "param.computation.method")
  # Delete attributes...
  for (i in seq_along(out)) {
    attr(out[[i]], "param.computation.method") <- NULL
  }

  do.call(returnParameters, c(
    list(out = out, distribution = distr),
    func = "parameters",
    lambdas = list(x$lambdas),
    trimmings = list(c(attr(x, "leftrim"), attr(x, "rightrim"))),
    tl.order = list(attr(x, "order")),
    param.computation.method = param.computation.method,
    attr(x, "source")
  ))
}

#' @method parameters.TLMoments data.frame
#' @export
parameters.TLMoments.data.frame <- function(x, distr, ...) {
  leftrim <- attr(x, "leftrim")
  rightrim <- attr(x, "rightrim")

  nam <- getFormulaSides(attr(x, "source")$formula)
  lambdas <- t(x$lambdas[!(names(x$lambdas) %in% nam$rhs)])
  out <- lapply(1L:ncol(lambdas), function(i) {
    param.est(lambdas[, i], leftrim, rightrim, distr, ...)
  })
  param.computation.method <- attr(out[[1]], "param.computation.method")
  out <- simplify2array(out)
  out <- cbind(x$lambdas[nam$rhs], as.data.frame(t(out)))

  do.call(returnParameters, c(
    list(out = out, distribution = distr),
    func = "parameters",
    lambdas = list(x$lambdas),
    trimmings = list(c(attr(x, "leftrim"), attr(x, "rightrim"))),
    tl.order = list(attr(x, "order")),
    param.computation.method = param.computation.method,
    attr(x, "source")
  ))
}

#' @export
print.parameters <- function(x, ...) {
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


param.est <- function(lambdas, leftrim, rightrim, distr, ...) {
  if (length(lambdas) < 2 && distr == "gpd" && "u" %in% names(list(...)))
    stop("Insufficient L-moments to calculate parameters. ")
  if (length(lambdas) < 3 && distr == "gpd" && !("u" %in% names(list(...))))
    stop("Insufficient L-moments to calculate parameters. ")
  if (length(lambdas) < 3 && distr %in% c("gev", "ln3"))
    stop("Insufficient L-moments to calculate parameters. ")
  if (length(lambdas) < 2 && distr == "gum")
    stop("Insufficient L-moments to calculate parameters. ")

  switch(
    distr,
    gev = gev.est(lambdas[1], lambdas[2], lambdas[3]/lambdas[2], leftrim, rightrim),
    gum = gum.est(lambdas[1], lambdas[2], leftrim, rightrim),
    gpd = gpd.est(lambdas[1], lambdas[2], lambdas[3]/lambdas[2], leftrim, rightrim, ...),
    ln3 = ln3.est(lambdas[1], lambdas[2], lambdas[3]/lambdas[2], leftrim, rightrim),
    stop("Not implemented. ")
  )
}

#### GUM ####
gum.est <- function(l1, l2, s, t) {
  switch(paste(s, t, sep = "-"),
         `0-0` = gum.tl00.est(l1, l2),
         `0-1` = gum.tl01.est(l1, l2),
         `1-0` = gum.tl10.est(l1, l2),
         `1-1` = gum.tl11.est(l1, l2),
         gum.numerical.est(l1, l2, s, t))
}
# TL(0,0)=L
gum.tl00.est <- function(l1, l2) {
  scale <- l2 / log(2)
  loc <- l1 - 0.577 * scale

  out <- setNames(c(loc, scale), c("loc", "scale"))
  attr(out, "param.computation.method") <- "analytical"
  return(out)
}
# TL(0,1), TL(1,0), TL(1,1) from Elamir, 2010: Optimal Choices for Trimming in Trimmed L-moment Method
gum.tl01.est <- function(l1, l2) {
  scale <- l2 / 0.431
  loc <- l1 + 0.116 * scale

  out <- setNames(c(loc, scale), c("loc", "scale"))
  attr(out, "param.computation.method") <- "analytical"
  return(out)
}
gum.tl10.est <- function(l1, l2) {
  scale <- l2 / 0.607
  loc <- l1 - 1.269 * scale

  out <- setNames(c(loc, scale), c("loc", "scale"))
  attr(out, "param.computation.method") <- "analytical"
  return(out)
}
gum.tl11.est <- function(l1, l2) {
  scale <- l2 / 0.353
  loc <- l1 - .459 * scale

  out <- setNames(c(loc, scale), c("loc", "scale"))
  attr(out, "param.computation.method") <- "analytical"
  return(out)
}
gum.numerical.est <- function(l1, l2, s, t) {
  #warning("Calculating numerical solution. ")

  # scale
  f1 <- function(scale) calcTLMom(2, s, t, qgum, loc = 0, scale = scale)[2] - l2
  i <- 0; while (sign(f1(10^-i)) == sign(f1(10^i))) i <- i+1
  u <- uniroot(f1, c(10^-i, 10^i))
  scale <- u$root
  # loc
  f2 <- function(loc) calcTLMom(2, s, t, qgum, loc = loc, scale = scale)[1] - l1
  i <- 0; while (sign(f2(-10^i)) == sign(f2(10^i))) i <- i+1
  u <- uniroot(f2, c(-10^i, 10^i))
  loc <- u$root

  out <- setNames(c(loc, scale), c("loc", "scale"))
  attr(out, "param.computation.method") <- "numerical"
  return(out)
}


#### GEV ####
gev.est <- function(l1, l2, t3, s, t) {
  switch(paste(s, t, sep = "-"),
         `0-0` = gev.tl00.est(l1, l2, t3),
         `0-1` = gev.tl01.est(l1, l2, t3),
         `0-2` = gev.tl02.est(l1, l2, t3),
         `1-1` = gev.tl11.est(l1, l2, t3),
         `1-0` = gev.tl10.est(l1, l2, t3),
         gev.numerical.est(l1, l2, t3, s, t))
}
# TL(0,0)=L der GEV
gev.tl00.est <- function(l1, l2, t3) {
  z <- 2/(3+t3) - log(2)/log(3)
  k <- 7.8590 * z + 2.9554 * z^2

  # Hosking (1990)
  gk <- gamma(1+k)
  al <- l2*k/((1-2^-k)*gk)
  xi <- l1 + al*(gk-1)/k

  out <- setNames(c(xi, al, -k), c("loc", "scale", "shape"))
  attr(out, "param.computation.method") <- "analytical"
  return(out)
}
# TL(0,1) der GEV
gev.tl01.est <- function(l1, l2, t3) {
  z <- 1/(2+t3) * 10/9 - (2*log(2)-log(3))/(3*log(3)-2*log(4))
  k <- 8.567394*z - 0.675969*z^2

  gk <- gamma(k)
  V3 <- (1/3)^k
  V2 <- (1/2)^k

  al <- l2 * 2/3 * 1/gk * 1/(V3 - 2*V2 + 1)
  xi <- l1 - al/k - al*gk*(V2 - 2)

  out <- setNames(c(xi, al, -k), c("loc", "scale", "shape"))
  attr(out, "param.computation.method") <- "analytical"
  return(out)
}
# TL(0,2) der GEV
gev.tl02.est <- function(l1, l2, t3) {

  z <- (t3+5/3) * 6/5 - (3*log(5)-8*log(4)+6*log(3))/(log(4)-3*log(3)+3*log(2))
  k <- -2.468959*z + 1.130074*z^2 - 0.635912*z^3

  gk <- gamma(k)
  V4 <- (1/4)^k
  V3 <- (1/3)^k
  V2 <- (1/2)^k

  al <- l2 * 2 * 1/gk * 1/(-4*V4 + 12*V3 - 12*V2 + 4)
  xi <- l1 - al/k - al*gk*(-V3 + 3*V2 - 3)

  out <- setNames(c(xi, al, -k), c("loc", "scale", "shape"))
  attr(out, "param.computation.method") <- "analytical"
  return(out)
}
# TL(1,0) der GEV
gev.tl10.est <- function(l1, l2, t3) {

  z <- (t3-40/3) * 9/20 - (log(3)-log(4))/(log(2)-log(3))
  k <- -9.066941*z - 3.374925*z^2 - 0.303208*z^3

  gk <- gamma(k)
  V3 <- (1/3)^k
  V2 <- (1/2)^k

  al <- l2 * 2 * 1/gk * 1/(3*V2 - 3*V3)
  xi <- l1 - al/k + al*gk*V2

  out <- setNames(c(xi, al, -k), c("loc", "scale", "shape"))
  attr(out, "param.computation.method") <- "analytical"
  return(out)
}
# TL(1,1) der GEV
gev.tl11.est <- function(l1, l2, t3) {

  z <- t3 * 9/20 + (log(3)-2*log(4)+log(5))/(log(2)-2*log(3)+log(4))
  k <- 25.3171*z - 91.5507*z^2 + 110.0626*z^3 - 46.5518*z^4

  gk <- gamma(k)
  V4 <- (1/4)^k
  V3 <- (1/3)^k
  V2 <- (1/2)^k

  al <- l2 * 1/gk * 1/(3*V2 - 6*V3 + 3*V4)
  xi <- l1 - al/k - al*gk*(-3*V2 + 2*V3)

  out <- setNames(c(xi, al, -k), c("loc", "scale", "shape"))
  attr(out, "param.computation.method") <- "analytical"
  return(out)
}
gev.numerical.est <- function(l1, l2, t3, s, t) {
  #warning("Calculating numerical solution. ")

  # shape
  f1 <- function(shape) {
    a <- calcTLMom(3, s, t, qgev, loc = 0, scale = 1, shape = shape)
    a[3]/a[2] - t3
  }
  i <- 0; while (sign(f1(-10^i)) == sign(f1(10^i))) i <- i+1
  u <- uniroot(f1, c(-10^i, 10^i))
  shape <- u$root
  # scale (positiv)
  f2 <- function(scale) calcTLMom(2, s, t, qgev, loc = 0, scale = scale, shape = shape)[2] - l2
  i <- 0; while (sign(f2(10^-i)) == sign(f2(10^i))) i <- i+1
  u <- uniroot(f2, c(10^-i, 10^i))
  scale <- u$root
  # loc
  f3 <- function(loc) calcTLMom(1, s, t, qgev, loc = loc, scale = scale, shape = shape)[1] - l1
  i <- 0; while (sign(f3(-10^i)) == sign(f3(10^i))) i <- i+1
  u <- uniroot(f3, c(-10^i, 10^i))
  loc <- u$root

  out <- setNames(c(loc, scale, shape), c("loc", "scale", "shape"))
  attr(out, "param.computation.method") <- "numerical"
  return(out)
}


#### LN3 ####
ln3.est <- function(l1, l2, t3, s, t) {
  switch(paste(s, t, sep = "-"),
         `0-0` = ln3.tl00.est(l1, l2, t3),
         ln3.numerical.est(l1, l2, t3, s, t))
}

ln3.tl00.est <- function(l1, l2, t3) {
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

  # 1. step, T3 <=> shape (positiv)
  f1 <- function(s) 6/sqrt(pi) * 1/erf(s/2) * integrate(function(x) erf(x/sqrt(3) * exp(-x^2)), 0, s/2)$val - t3
  i <- 1; while (sign(f1(10^-i)) == sign(f1(sqrt(exp(i))))) i <- i+1
  u <- uniroot(f1, c(10^-i, sqrt(exp(i))), tol = .Machine$double.eps^0.5, check.conv = TRUE)
  shape <- u$root

  # 2. scale and loc
  scale <- log(l2 / erf(shape/2)) - shape^2/2
  loc <- l1 - exp(scale + shape^2/2)

  out <- setNames(c(loc, scale, shape), c("loc", "scale", "shape"))
  attr(out, "param.computation.method") <- "analytical"
  return(out)
}

ln3.numerical.est <- function(l1, l2, t3, s, t) {
  #warning("Calculating numerical solution. ")

  # 1. step, T3 <=> shape (positiv)
  f1 <- function(shape) {
    a <- calcTLMom(3, s, t, qln3, loc = 0, scale = 1, shape = shape)
    a[3]/a[2] - t3
  }
  i <- 1; while (sign(f1(10^-i)) == sign(f1(sqrt(exp(i))))) i <- i+1
  u <- uniroot(f1, c(10^-i, sqrt(exp(i))))
  shape <- u$root
  # 2. step, L2 <=> scale
  f2 <- function(scale) calcTLMom(2, s, t, qln3, loc = 0, scale = scale, shape = shape)[2] - l2
  i <- 0; while (sign(f2(-10^i)) == sign(f2(10^i))) i <- i+1
  u <- uniroot(f2, c(-10^i, 10^i))
  scale <- u$root
  # 3. step L1 <=> loc
  f3 <- function(loc) calcTLMom(1, s, t, qln3, loc = loc, scale = scale, shape = shape)[1] - l1
  i <- 0; while (sign(f3(-10^i)) == sign(f3(10^i))) i <- i+1
  u <- uniroot(f3, c(-10^i, 10^i))
  loc <- u$root

  out <- setNames(c(loc, scale, shape), c("loc", "scale", "shape"))
  attr(out, "param.computation.method") <- "numerical"
  return(out)
}


#### GPD ####
gpd.est <- function(l1, l2, t3, s, t, ...) {
  switch(paste(s, t, sep = "-"),
         `0-0` = gpd.tl00.est(l1, l2, t3, ...),
         `0-1` = gpd.tl01.est(l1, l2, t3, ...),
         `0-2` = gpd.tl02.est(l1, l2, t3, ...),
         gpd.numerical.est(l1, l2, t3, s, t, ...))
}
# GPD
gpd.tl00.est <- function(l1, l2, t3, u = NULL) {

  if (!is.null(u) & is.numeric(u)) {
    shape <- -(l1-u)/l2 + 2
    scale <- (l1-u) * (-shape+1)
    loc <- u
  } else {
    shape <- (3*t3 - 1)/(t3 + 1)
    scale <- l2 * (shape-1) * (shape-2)
    loc <- l1 + scale/(shape-1)
  }

  out <- setNames(c(loc, scale, shape), c("loc", "scale", "shape"))
  attr(out, "param.computation.method") <- "analytical"
  return(out)
}

gpd.tl01.est <- function(l1, l2, t3, u = NULL) {

  if (!is.null(u) & is.numeric(u)) {
    # from Hosking, 2007.
    shape <- -3/2 * (l1-u)/l2 + 3
    scale <- (l1-u) * (-shape+2)
    loc <- u
  } else {
    # from Fischer, 2015.
    shape <- (36*t3 - 8)/(9*t3 + 8)
    scale <- 2/3 * l2 * (shape-2) * (shape-3)
    loc <- l1 + scale/(shape-2)
  }

  out <- setNames(c(loc, scale, shape), c("loc", "scale", "shape"))
  attr(out, "param.computation.method") <- "analytical"
  return(out)
}

gpd.tl02.est <- function(l1, l2, t3, u = NULL) {

  if (!is.null(u) & is.numeric(u)) {
    # from Hosking, 2007.
    shape <- -2 * (l1-u)/l2 + 4
    scale <- (l1-u) * (-shape+3)
    loc <- u
  } else {
    # from Fischer, 2015.
    shape <- (30*t3 - 5)/(6*t3 + 5)
    scale <- 1/2 * l2 * (shape-3) * (shape-4)
    loc <- l1 + scale/(shape-3)
  }

  out <- setNames(c(loc, scale, shape), c("loc", "scale", "shape"))
  attr(out, "param.computation.method") <- "analytical"
  return(out)
}

gpd.numerical.est <- function(l1, l2, t3, s, t, u = NULL) {
  #warning("Calculating numerical solution. ")

  if (!is.null(u) & is.numeric(u)) {

    # shape
    f1u <- function(shape, u) {
      a <- calcTLMom(2, s, t, qgpd, loc = u, scale = 1, shape = shape)
      a[2] - l2
    }
    i <- 0; while (sign(f1u(-10^i, u)) == sign(f1u(10^i, u))) i <- i+1
    ur <- uniroot(f1u, c(-10^i, 10^i), u = u)
    shape <- ur$root
    # scale
    f2u <- function(scale, u) calcTLMom(1, s, t, qgpd, loc = u, scale = scale, shape = shape)[1] - l1
    i <- 0; while (sign(f2u(10^-i, u)) == sign(f2u(10^i, u))) i <- i+1
    ur <- uniroot(f2u, c(10^-i, 10^i), u = u)
    scale <- ur$root
    # loc
    loc <- u

  } else {

    # shape
    f1 <- function(shape) {
      a <- calcTLMom(3, s, t, qgpd, loc = 0, scale = 1, shape = shape)
      a[3]/a[2] - t3
    }
    i <- 0; while (sign(f1(-10^i)) == sign(f1(10^i))) i <- i+1
    ur <- uniroot(f1, c(-10^i, 10^i))
    shape <- ur$root
    # scale
    f2 <- function(scale) calcTLMom(2, s, t, qgpd, loc = 0, scale = scale, shape = shape)[2] - l2
    i <- 0; while (sign(f2(10^-i)) == sign(f2(10^i))) i <- i+1
    ur <- uniroot(f2, c(10^-i, 10^i))
    scale <- ur$root
    # loc
    f3 <- function(loc) calcTLMom(1, s, t, qgpd, loc = loc, scale = scale, shape = shape)[1] - l1
    i <- 0; while (sign(f3(-10^i)) == sign(f3(10^i))) i <- i+1
    ur <- uniroot(f3, c(-10^i, 10^i))
    loc <- ur$root

  }

  out <- setNames(c(loc, scale, shape), c("loc", "scale", "shape"))
  attr(out, "param.computation.method") <- "numerical"
  return(out)
}
