#' @title
#' Estimate the covariance matrix of TL-moments estimations
#' @description
#' Internal function. Use \link{est_cov}. Description not done yet.
#'
#' @param x numeric vector or matrix containing data OR an object of TLMoments.
#' @param leftrim,rightrim integer indicating lower and upper trimming parameters, have to be non-negative integers.
#' @param order numeric vector giving the orders that are returned (default is first three L-moments).
#' @param distr character of length 1 giving the distribution if parametric assumption should be used.
#' @param set.n hypothetical data length n if theoretical values are given.
#' @param lambda.cov boolean, if TRUE (default) TL-moment estimation covariance matrix is calculated.
#' @param ratio.cov boolean, if TRUE (default) TL-moment-ratio estimation covariance matrix is calculated.
#' @param reg.weights numeric vector of weights for regionalized TLMoments.
#' @param ... additional arguments.
#'
#' @return a list of numeric matrices (if \code{lambda.cov} and \code{ratio.cov} are TRUE (default)), or a single matrix.
#'
#' @examples
#' ### Numeric vectors
#' x <- rgev(500, loc = 10, scale = 5, shape = .1)
#'
#' est_tlmcov(x)
#' est_tlmcov(x, order = 2:3)
#' est_tlmcov(x, rightrim = 1, order = 4:5)
#' # cov(t(replicate(10000,
#' #   TLMoments(rgev(500, loc = 10, scale = 5, shape = .1))$lambdas)
#' # ))
#' # cov(t(replicate(10000,
#' #   TLMoments(rgev(500, loc = 10, scale = 5, shape = .1))$ratios)
#' # ))
#'
#' est_tlmcov(x, ratio.cov = FALSE)
#' est_tlmcov(x, lambda.cov = FALSE)
#'
#' est_tlmcov(x, distr = "gev")
#'
#' est_tlmcov(x, leftrim = 0, rightrim = 1)
#' # cov(t(replicate(10000,
#' #  TLMoments(rgev(500, loc = 10, scale = 5, shape = .1), 0, 1, 3)$lambdas
#' # )))
#' # cov(t(replicate(10000,
#' #  TLMoments(rgev(500, loc = 10, scale = 5, shape = .1), 0, 1, 3)$ratios
#' # )))
#'
#' ### Numeric matrices
#' x <- matrix(rgev(600), nc = 3)
#'
#' est_tlmcov(x)
#' est_tlmcov(x, order = 3:4)
#' # cov(t(replicate(10000,
#' #   as.vector(TLMoments(matrix(rgev(600), nc = 3))$lambdas[3:4, ])
#' # )))
#' # cov(t(replicate(10000,
#' #   as.vector(TLMoments(matrix(rgev(600), nc = 3))$ratios[3:4, ])
#' # )))
#'
#' est_tlmcov(x, ratio.cov = FALSE)
#' est_tlmcov(x, lambda.cov = FALSE)
#'
#' TLMoments:::est_tlmcov(x, order = 2:3, distr = "gev")
#' # cov(t(replicate(10000,
#' #   as.vector(TLMoments(matrix(rgev(600), nc = 3))$lambdas[2:3, ])
#' # )))
#' # cov(t(replicate(10000,
#' #   as.vector(TLMoments(matrix(rgev(600), nc = 3))$ratios[2:3, ])
#' # )))
#'
#' ### TLMoments-object (theoretical calculation)
#' tlm <- TLMoments(as.parameters(loc = 10, scale = 5, shape = .1, distr = "gev"), 0, 1)
#' est_tlmcov(tlm, distr = "gev", set.n = 100)
#' est_tlmcov(tlm, distr = "gev", set.n = 100, ratio.cov = FALSE)
#' est_tlmcov(tlm, distr = "gev", set.n = 100, lambda.cov = FALSE)
#'
#' @rdname est_tlmcov
#' @export
est_tlmcov <- function(x,
                       leftrim = 0L,
                       rightrim = 0L,
                       order = 1:3,
                       distr = NULL,
                       lambda.cov = TRUE,
                       ratio.cov = TRUE,
                       ...) {
  # TODO: General error catches
  if (!are.integer.like(leftrim, rightrim))
    stop("leftrim and rightrim have to be integer types!")
  if (!is.logical(lambda.cov) || !is.logical(ratio.cov))
    stop("Arguments \"lambda.cov\" and \"ratio.cov\" must be logical/boolean.")
  if (!is.null(distr) && distr != "gev")
    stop("distr has to be NULL (nonparametric) or \"gev\"")

  UseMethod("est_tlmcov")
}

#' @rdname est_tlmcov
#' @method est_tlmcov numeric
#' @export
est_tlmcov.numeric <- function(x,
                               leftrim = 0L,
                               rightrim = 0L,
                               order = 1:3,
                               distr = NULL,
                               lambda.cov = TRUE,
                               ratio.cov = TRUE,
                               ...) {

  maxk <- (max(order)+leftrim+rightrim)
  pwmcov <- est_pwmcov(x, 0:(maxk-1), distr = distr, distr.trim = c(leftrim, rightrim))

  Z <- sapply(0:(maxk-1), function(k) sapply(0:(max(order)-1), function(r) {
    z_C(r, k, leftrim, rightrim)
  }))
  lambdacov <- (Z %*% pwmcov %*% t(Z))
  rownames(lambdacov) <- colnames(lambdacov) <- buildNames("L", 1L:max(order))

  if (max(order) >= 2 && ratio.cov) {
    A <- CovLambdaToTau(TLMoments(x, leftrim = leftrim, rightrim = rightrim, max.order = max(order), na.rm = TRUE)$lambdas)
    taucov <- t(A) %*% lambdacov %*% A
    rownames(taucov) <- colnames(taucov) <- buildNames("T", 2L:max(order))
  } else {
    taucov <- NULL
  }

  if (lambda.cov && ratio.cov) {
    out <- list(
      lambdas = lambdacov[paste0("L", order), paste0("L", order), drop = FALSE],
      ratios = taucov[paste0("T", order[order!=1]), paste0("T", order[order!=1]), drop = FALSE]
    )
  } else if (lambda.cov) {
    out <- lambdacov[paste0("L", order), paste0("L", order), drop = FALSE]
  } else if (ratio.cov) {
    out <- taucov[paste0("T", order[order!=1]), paste0("T", order[order!=1]), drop = FALSE]
  } else stop("Invalid arguments given. ")

  attr(out, "distribution") <- distr
  attr(out, "trimmings") <- c(leftrim, rightrim)
  attr(out, "n") <- length(x)
  out
}

#' @rdname est_tlmcov
#' @method est_tlmcov matrix
#' @export
est_tlmcov.matrix <- function(x,
                              leftrim = 0L,
                              rightrim = 0L,
                              order = 1:3,
                              distr = NULL,
                              lambda.cov = TRUE,
                              ratio.cov = TRUE,
                              reg.weights = NULL,
                              ...) {

  if (!is.null(reg.weights) && length(reg.weights) != ncol(x))
    stop("Invalid reg.weights. ")

  maxk <- (max(order)+leftrim+rightrim)
  pwmcov <- est_pwmcov(x, 0:(maxk-1), distr = distr, distr.trim = c(leftrim, rightrim))

  z <- sapply(0:(maxk-1), function(k) sapply(0:(max(order)-1), function(r) {
    z_C(r, k, leftrim, rightrim)
  }))

  Z <- blockdiag(z, ncol(x))
  lambdacov <- (Z %*% pwmcov %*% t(Z))
  lnames <- buildNames("L", 1L:max(order), 1L:ncol(x))
  lidx <- grep(paste0("L[", paste0(order, collapse = "|"), "]"), x = lnames)
  rownames(lambdacov) <- colnames(lambdacov) <- lnames

  if (max(order) >= 2 && is.null(reg.weights)) {
    A <- lapply(1:ncol(x), function(i) CovLambdaToTau(TLMoments(x[, i], leftrim = leftrim, rightrim = rightrim, max.order = max(order), na.rm = TRUE)$lambdas))
    A <- blockdiag_list(A)
    taucov <- t(A) %*% lambdacov %*% A
    tnames <- buildNames("T", 2L:max(order), 1L:ncol(x))
    tidx <- grep(paste0("T[", paste0(order, collapse = "|"), "]"), x = tnames)
    rownames(taucov) <- colnames(taucov) <- tnames
  }

  # regionalized:
  if (!is.null(reg.weights)) {
    R <- do.call(cbind, lapply(reg.weights, function(x) diag(rep(x, length(order)))))
    lambdacov <- R %*% lambdacov %*% t(R)
    lnames <- buildNames("L", 1L:max(order))
    lidx <- grep(paste0("L[", paste0(order, collapse = "|"), "]"), x = lnames)
    rownames(lambdacov) <- colnames(lambdacov) <- lnames

    if (max(order) >= 2) {
      A <- CovLambdaToTau(regionalize(TLMoments(x, leftrim = leftrim, rightrim = rightrim, max.order = max(order), na.rm = TRUE), reg.weights)$lambdas)
      taucov <- t(A) %*% lambdacov %*% A
      tnames <- buildNames("T", 2L:max(order))
      tidx <- grep(paste0("T[", paste0(order, collapse = "|"), "]"), x = tnames)
      rownames(taucov) <- colnames(taucov) <- tnames
    }
  }

  if (lambda.cov && ratio.cov) {
    out <- list(
      lambdas = lambdacov[lidx, lidx],
      ratios = taucov[tidx, tidx]
    )
  } else if (lambda.cov) {
    out <- lambdacov[lidx, lidx]
  } else if (ratio.cov) {
    out <- taucov[tidx, tidx]
  } else stop("Invalid arguments given. ")

  attr(out, "distribution") <- distr
  attr(out, "trimmings") <- c(leftrim, rightrim)
  attr(out, "n") <- apply(x, 2, function(y) sum(!is.na(y)))
  out
}

#' @rdname est_tlmcov
#' @method est_tlmcov TLMoments
#' @export
est_tlmcov.TLMoments <- function(x,
                                 leftrim = attr(x, "leftrim"),
                                 rightrim = attr(x, "rightrim"),
                                 order = attr(x, "order"),
                                 distr = NULL,
                                 lambda.cov = TRUE,
                                 ratio.cov = TRUE,
                                 set.n = NA,
                                 ...) {

  if (!inherits(x$lambdas, "numeric"))
    stop("x must be a numeric TLMoments-object. ")

  if (!any(attr(x, "source")$func %in% c("as.PWMs", "as.TLMoments", "as.parameters")))
    stop("est_tlmcov.TLMoments only for theoretical values. ")

  if (is.null(distr))
    distr <- attr(attr(x, "source")$data, "distribution")
    if (is.null(distr)) stop("distr argument must be given.")

  if (is.na(set.n) | !is.numeric(set.n)) {
    #warning("Missing or invalid set.n argument. Giving results for n = 100. ")
    n <- 100
  } else { n <- set.n }

  maxk <- (max(order)+leftrim+rightrim)
  p <- parameters(x, distr)
  r <- parametricPWMCov(distr, 0:(maxk-1), scale = p["scale"], shape = p["shape"])
  pwmcov <- r / n

  Z <- Z_C(max(order), leftrim, rightrim)
  lambdacov <- (Z %*% pwmcov %*% t(Z))
  rownames(lambdacov) <- colnames(lambdacov) <- paste0("L", 1:max(order))

  if (max(order) >= 2 && ratio.cov) {
    A <- CovLambdaToTau(x$lambdas[1:max(order)])
    taucov <- t(A) %*% lambdacov %*% A
    rownames(taucov) <- colnames(taucov) <- paste0("T", 2:max(order))
  }

  if (lambda.cov && ratio.cov) {
    out <- list(
      lambdas = lambdacov[paste0("L", order), paste0("L", order)],
      ratios = taucov[paste0("T", order[order!=1]), paste0("T", order[order!=1])]
    )
  } else if (lambda.cov) {
    out <- lambdacov[paste0("L", order), paste0("L", order)]
  } else if (ratio.cov) {
    out <- taucov[paste0("T", order[order!=1]), paste0("T", order[order!=1])]
  } else stop("Invalid arguments given. ")

  attr(out, "distribution") <- distr
  attr(out, "trimmings") <- c(leftrim, rightrim)
  attr(out, "n") <- n
  out
}
