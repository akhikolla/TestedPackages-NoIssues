read.dsm.matrix <- function (file, format=c("word2vec"), encoding="UTF-8", batchsize=1e6, verbose=FALSE) {
  format <- match.arg(format)
  
  if (is.character(file)) {
    file <- file(file, "r", encoding=encoding)
    on.exit(close(file))
  }
  else {
    if (!inherits(file, "connection")) stop("'file' must be a character string or connection")
    if (!isOpen(file, "r")) {
      open(file, "r")
      on.exit(close(file))
    }
  }
  
  if (format == "word2vec") {
    info <- read.table(file, header=FALSE, quote="", comment.char="", colClasses="integer", col.names=c("nR", "nC"), nrows=1)
    if (!(info$nR > 0 && info$nC > 0)) stop("format error in 'word2vec' file")
    if (as.double(info$nR) * info$nC > 2147483647) stop(sprintf("%d x %d matrix is too large for R", info$nR, info$nC))
    batch.lines <- max(round(batchsize / info$nC), 10) # how many matrix rows to read per batch
    col.classes <- c("character", rep("double", info$nC))
    M <- matrix(0, info$nR, info$nC) # pre-allocate the data matrix
    row.names <- character(info$nR)
    n.total <- 0
    if (verbose) {
      pb <- txtProgressBar(0, info$nR, style=3)
      on.exit(close(pb), add=TRUE)
    }
    while (n.total < info$nR) {
      A <- as.matrix(read.table(file, header=FALSE, quote="", comment.char="", na.strings=c(), colClasses=col.classes, row.names=1, nrows=batch.lines))
      if (ncol(A) != info$nC) stop(sprintf("format error in 'word2vec' file -- expected %d-dim vector, but got %d columns", info$nC, ncol(A)))
      nR <- nrow(A)
      if (nR < 1) stop(sprintf("read error in 'word2vec' file -- expecting to read %d further rows", info$nR - n.total))
      if (n.total + nR > info$nR) stop("format error in 'word2vec' file -- too many rows")
      idx <- seq(n.total + 1, n.total + nR)
      M[idx, ] <- A
      row.names[idx] <- rownames(A)
      n.total <- n.total + nR
      if (verbose) setTxtProgressBar(pb, n.total)
    }
    structure(M, dimnames=list(row.names, NULL))
  }
  else stop(sprintf("internal error -- unsupported format '%s'", format))
}
