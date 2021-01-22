t.dsm <- function (x) {
  info <- check.dsm(x)
  tx <- list(rows=x$cols, cols=x$rows, globals=x$globals)
  if (info$M$ok) tx$M <- t(x$M)
  if (info$S$ok) tx$S <- t(x$S)
  class(tx) <- c("dsm", "list")
  return(tx)
}
