read.dsm.triplet <- function (filename, freq=FALSE, value.first=FALSE, tokens=FALSE,
                              rowinfo=NULL, rowinfo.header=NULL, colinfo=NULL, colinfo.header=NULL, N=NA, span.size=1,
                              sep="\t", quote="", nmax=-1, sort=FALSE, encoding=getOption("encoding"), verbose=FALSE) {
  if (verbose) cat(sprintf("Loading DSM triplets from '%s' ... ", filename))
  is.pipe <- grepl("\\|\\s*$", filename, perl=TRUE)
  if (is.pipe) {
    filename <- sub("\\s*\\|\\s*$", "", filename, perl=TRUE)
    fh <- pipe(filename) # don't set encoding= parameter because .read.delim.raw will read binary data and convert later
  } else {
    fh <- file(filename) # don't set encoding= parameter
  }
  
  if (tokens) {
    if (freq || value.first) warning("freq= and value.first= options are ignored with tokens=TRUE")
    triplets <- .read.delim.raw(fh, header=FALSE, colClasses=c(l1="character", l2="character"), sep=sep, quote=quote, nrows=nmax, encoding=encoding)
    ## I prefer iotools::read.delim.raw over readr::read_delim because:
    ##  - readr has many "expensive" dependencies (esp. Boost in package 'BH')
    ##  - readr doesn't support the default "native.enc" encoding and cannot read from all types of connections
    ##  - iotools is slightly faster and leaner (less memory overhead) than readr
    ##  - unfortunately, read.delim.raw doesn't convert character encodings at all, but we work around this with our own .read.delim.raw
    ## Alternative version using readr::read_delim:
    ##   triplets <- read_delim(fh, sep, col_names=c("l1", "l2"), col_types="cc", locale=locale(encoding=encoding), n_max=nmax, quote=quote, comment="", na=character(), escape_double=FALSE, escape_backslash=FALSE, progress=FALSE)
    triplets$val <- 1
    freq <- TRUE
  } else {
    col.types <- if (value.first) c(val="numeric", l1="character", l2="character") else c(l1="character", l2="character", val="numeric")
    triplets <- .read.delim.raw(fh, header=FALSE, colClasses=col.types, sep=sep, quote=quote, nrows=nmax, encoding=encoding)
    ## Alternative version using readr::read_delim:
    ##   col.types <- if (value.first) "dcc" else "ccd"
    ##   col.names <- if (value.first) c("val", "l1", "l2") else c("l1", "l2", "val")
    ##   triplets <- read_delim(fh, sep, col_names=col.names, col_types=col.types, locale=locale(encoding=encoding), n_max=nmax, quote=quote, comment="", na=character(), escape_double=FALSE, escape_backslash=FALSE, progress=FALSE)
  }
  ## close(fh) not needed because .read.delim.raw automatically opens and closes the connection
  if (verbose) cat(sprintf("%.2fM %s\n", length(triplets$l1) / 1e6, if (tokens) "tokens" else "items"))

  ## read external marginal frequencies or other information on targets and features
  have.rowinfo <- !is.null(rowinfo)
  have.colinfo <- !is.null(colinfo)
  have.N <- !is.na(N)
  if ((have.rowinfo || have.colinfo) && freq) {
    if (!(have.rowinfo && have.colinfo && have.N)) stop("need rowinfo=, colinfo= and N= for external marginal frequencies")
  }
  
  if (have.rowinfo) {
    if (verbose) cat(sprintf(" - loading target information from '%s'\n", rowinfo))
    if (is.null(rowinfo.header)) rowinfo.header <- TRUE
    rowinfo.tbl <- .read.delim.raw(file(rowinfo), header=rowinfo.header, sep=sep, quote=quote, encoding=encoding)
    if (!("term" %in% colnames(rowinfo.tbl))) stop("rowinfo= must specify feature types in column 'term'")
    have.f <- "f" %in% colnames(rowinfo.tbl)
    if (freq && !have.f) stop("rowinfo= must include marginal frequencies (column 'f') if freq=TRUE or tokens=TRUE")
    if (have.f && !missing(span.size)) rowinfo.tbl$f <- span.size * rowinfo.tbl$f # adjust row marginals for span size
    ## dsm() constructor below will check that all target terms are included in the table and add nnzero counts if necessary
  } 
  else rowinfo.tbl <- NULL

  if (have.colinfo) {
    if (verbose) cat(sprintf(" - loading feature information from '%s'\n", colinfo))
    if (is.null(colinfo.header)) colinfo.header <- TRUE
    colinfo.tbl <- .read.delim.raw(file(colinfo), header=colinfo.header, sep=sep, quote=quote, encoding=encoding)
    if (!("term" %in% colnames(colinfo.tbl))) stop("colinfo= must specify target types in column 'term'")
    have.f <- "f" %in% colnames(colinfo.tbl)
    if (freq && !have.f) stop("colinfo= must include marginal frequencies (column 'f') if freq=TRUE or tokens=TRUE")
  }
  else colinfo.tbl <- NULL
  
  dsm(target=triplets$l1, feature=triplets$l2, score=triplets$val, raw.freq=freq, rowinfo=rowinfo.tbl, colinfo=colinfo.tbl, N=N, sort=sort, verbose=verbose)
}
