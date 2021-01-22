print.dsm <- function (x, ...) {
  info <- check.dsm(x)
  cat(sprintf("Distributional Semantic Model with %d rows x %d columns\n", info$nrow, info$ncol))
  n.cells <- prod(info$nrow, info$ncol)
  qformat <- function (y) {
    if (is.na(y)) return("N/A")
    if (y > 9.9999e9) return(sprintf("%.1fG", y/1e9))
    if (y > 9.9999e6) return(sprintf("%.1fM", y/1e6))
    if (y > 50e3) return(sprintf("%.1fk", y/1e3))
    return(sprintf("%d", y))
  }
  if (info$M$ok) {
    cat("* raw co-occurrence matrix M available\n")
    if (info$M$sparse) {
      n.nz <- signcount(x$M, "nnzero")
      cat(sprintf("  - sparse matrix with %s / %s nonzero entries (fill rate = %.2f%%)\n", qformat(n.nz), qformat(n.cells), 100 * n.nz / n.cells))
    } else {
      cat(sprintf("  - dense matrix with %s cells\n", qformat(n.cells)))
    }
    if (info$M$canonical) cat("  - in canonical format\n")
    if (isTRUE(info$M$nonneg)) cat("  - known to be non-negative\n")
    if (!is.na(info$N)) cat(sprintf("  - sample size of underlying corpus: %s tokens\n", qformat(info$N)))
  }
  if (info$S$ok) {
    cat("* scored matrix S available\n")
    if (info$S$sparse) {
      n.nz <- signcount(x$S, "nnzero")
      cat(sprintf("  - sparse matrix with %s / %s nonzero entries (fill rate = %.2f%%)\n", qformat(n.nz), qformat(n.cells), 100 * n.nz / n.cells))
    } else {
      cat(sprintf("  - dense matrix with %s cells\n", qformat(n.cells)))
    }
    if (info$S$canonical) cat("  - in canonical format\n")
    if (isTRUE(info$S$nonneg)) cat("  - known to be non-negative\n")
  }
  invisible(x)
}
