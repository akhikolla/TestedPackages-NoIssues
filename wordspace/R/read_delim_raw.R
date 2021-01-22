## replacement for iotools::read.delim.raw, which re-encodes input file after reading
## - most options are fixed and column names and classes have to be specified by the user
## - fh must be a connection that can be properly opened in binary mode (with automatic decompression if necessary)
## - column names can be read from file (header=TRUE), specified as character vector or implicitly as names of colClasses
## - colClasses must be given as character strings; if unspecified, they are inferred from the first lines of the file
## - unless encoding is "" or "native.enc", input will be converted to UTF-8 and marked as such
.read.delim.raw <- function (fh, header=TRUE, colClasses=NULL, sep="\t", quote="", nrows=-1L, nrowsClasses=25L, encoding=getOption("encoding")) {
  is.native <- encoding == "" || encoding == "native.enc"
  is.utf8 <- any(grepl("^utf-?8$", encoding, ignore.case=TRUE))
  with.header <- isTRUE(header)
  
  ## read entire file as raw vector (may want to make block size n= larger than default of 64k bytes)
  bindata <- readAsRaw(fh)
  if (!is.native && !is.utf8) {
    ## use iconv to recode the raw vector into UTF-8
    bindata <- iconv(list(bindata), from=encoding, to="UTF-8", toRaw=TRUE)[[1]]
  }
  
  ## if colClasses aren't specified, try to infer them from the first nrowsClasses lines
  if (is.null(colClasses)) {
    subset <- mstrsplit(bindata, sep=sep, quote=quote, nsep=NA, strict=TRUE, nrows=nrowsClasses, skip=with.header)
    colClasses <- sapply(1:ncol(subset), function (i) {
      type <- class(type.convert(subset[, i], as.is=TRUE))
      if (type == "integer") "numeric" else type # integer marginal frequencies can be problematic
    })
  }
  
  ## read header line
  if (with.header) {
    header <- mstrsplit(bindata, sep=sep, quote=quote, nsep=NA, strict=TRUE, nrows=1L)
  } 
  
  ## add header names to colClasses
  if (is.character(header)) {
    if (length(header) != length(colClasses)) stop("numer of items in header doesn't match number of table columns")
    names(colClasses) <- header
  }
  else {
    if (is.null(names(colClasses))) names(colClasses) <- sprintf("V%d", seq_along(colClasses))
  }
  
  ## now parse the raw table data into a data frame
  res <- dstrsplit(bindata, colClasses, sep=sep, quote=quote, nsep=NA, strict=TRUE, nrows=nrows, skip=with.header)
  if (!is.native) {
    ## mark all character variables as UTF-8 (unless read with native encoding)
    for (i in seq_along(colClasses)) {
      if (colClasses[i] == "character") Encoding(res[[i]]) <- "UTF-8"
    }
  }
  res
}
