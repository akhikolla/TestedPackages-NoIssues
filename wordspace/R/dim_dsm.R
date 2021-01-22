dim.dsm <- function (x) {
  res <- check.dsm(x)
  c(res$nrow, res$ncol)
}