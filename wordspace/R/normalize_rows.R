## returns scaling exponent required for normalization (usually 1)
.check.normalize <- function (method, p) {
  if (method == "minkowski") {
    if (p == 0) stop("Hamming length (p = 0) cannot be normalized")
    if (p < .05) stop("reliable normalization not possible for Minkowski norm with small p (< .05)")
    if (p < 1) return(1 / p)
  }
  return(1)
}

normalize.rows <- function (M, method = "euclidean", p = 2, ..., tol = 1e-6, inplace = FALSE) {
  scale <- .check.normalize(method, p)
  norms <- rowNorms(M, method=method, p=p, ...)
  lambda <- ifelse(norms < tol, 0, 1 / norms)  # any rows with norm < tol are set to 0
  if (scale == 1) scaleMargins(M, rows = lambda, duplicate=!inplace) else scaleMargins(M, rows = lambda ^ scale, duplicate=!inplace)
}

normalize.cols <- function (M, method = "euclidean", p = 2, ..., tol = 1e-6, inplace = FALSE) {
  scale <- .check.normalize(method, p)
  norms <- colNorms(M, method=method, p=p, ...)
  lambda <- ifelse(norms < tol, 0, 1 / norms)  # any cols with norm < tol are set to 0
  if (scale == 1) scaleMargins(M, cols = lambda, duplicate=!inplace) else scaleMargins(M, cols = lambda ^ scale, duplicate=!inplace)
}
