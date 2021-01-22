as.dsm.TermDocumentMatrix <- function (obj, ..., verbose=FALSE) {
  if (! inherits(obj, c("TermDocumentMatrix", "DocumentTermMatrix"))) stop("argument must be a 'tm' object of class 'TermDocumentMatrix' or 'DocumentTermMatrix'")
  if (any(is.na(obj$v))) stop("missing values in 'tm' matrix are not supported by wordspace package")
  M <- sparseMatrix(i=obj$i, j=obj$j, x=obj$v, dims=c(obj$nrow, obj$ncol), dimnames=obj$dimnames)
  raw.freq <- "tf" %in% attr(obj, "Weighting")
  dsm(M, raw.freq=raw.freq, sort=FALSE, verbose=verbose)
}

as.dsm.DocumentTermMatrix <- as.dsm.TermDocumentMatrix
