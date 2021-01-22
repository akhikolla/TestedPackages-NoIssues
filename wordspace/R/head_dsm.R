head.dsm <- function (x, n=6L, k=n, ...) {
  info <- check.dsm(x)
  k <- min(k, info$ncol) # do this first, so change in n doesn't affect the default
  n <- min(n, info$nrow)
  M <- if (info$S$ok) x$S else x$M
  M[1:n, 1:k]
}
