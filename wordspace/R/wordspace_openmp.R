wordspace.openmp <- function (threads=NULL) {
  if (!is.null(threads)) {
    if (!(is.numeric(threads) && length(threads) == 1 && threads >= 1)) stop("argument threads= must be a single integer >= 1")
    CPP_set_openmp_threads(threads)
  }
  res <- CPP_get_openmp_threads()
  if (is.null(threads)) res else invisible(res)
}
