head.dist.matrix <- function (x, n=6L, k=n, ...) {
  n <- min(n, nrow(x))
  k <- min(k, ncol(x))
  x[1:n, 1:k]
}

