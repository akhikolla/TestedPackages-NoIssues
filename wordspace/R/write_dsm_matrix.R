write.dsm.matrix <- function (x, file, format=c("word2vec"), round=FALSE, encoding="UTF-8", batchsize=1e6, verbose=FALSE) {
  format <- match.arg(format)
  M <- find.canonical.matrix(x)
  info <- dsm.is.canonical(M, nonneg.check=FALSE)
  nR <- nrow(M)
  nC <- ncol(M)

  if (format == "word2vec") {
    if (info$sparse) stop("'word2vec' format cannot be used for a sparse matrix")
    idx <- grepl("\\s", rownames(M), perl=TRUE)
    if (any(idx)) stop(sprintf("row labels must not contain whitespace in word2vec format (%d problematic labels found)", sum(idx)))
    idx <- duplicated(rownames(M))
    if (any(idx)) stop(sprintf("duplicate row labels are not allowed in word2vec format (%d problematic labels found)", sum(idx)))
  }
  
  if (is.character(file)) {
    file <- file(file, "w", encoding=encoding)
    on.exit(close(file))
  }
  else {
    if (!inherits(file, "connection")) stop("'file' must be a character string or connection")
    if (!isOpen(file, "w")) {
      open(file, "w")
      on.exit(close(file))
    }
  }
  
  if (format == "word2vec") {
    batch.lines <- max(round(batchsize / nC), 10) # how many matrix rows to write per batch
    if (verbose) {
      pb <- txtProgressBar(0, nR, style=3)
      on.exit(close(pb), add=TRUE)
    }
    cat(sprintf("%d %d\n", nrow(M), ncol(M)), file=file)
    for (i in seq(1, nR, batch.lines)) {
      k <- min(i + batch.lines - 1, nR)
      M.batch <- M[i:k, , drop=FALSE]
      if (!isFALSE(round)) M.batch <- round(M.batch, round)
      write.table(M.batch, file=file, quote=FALSE, sep=" ", row.names=TRUE, col.names=FALSE)
      if (verbose) setTxtProgressBar(pb, k)
    }
  }
  else stop(sprintf("internal error -- unsupported format '%s'", format))
}
