## Test reading data sets in different file formats, based on the example DSM_TermContext matrix
context("Input file formats")
library(wordspace)

## example data: DSM_TermContext
TC <- DSM_TermContext        # order rows and columns alphabetically for comparison with input files
idxR <- order(TC$rows$term)
idxC <- order(TC$cols$term)
TC$M <- TC$M[idxR, idxC]
TC$rows <- TC$rows[idxR, ]
TC$cols <- TC$cols[idxC, ]
invisible(check.dsm(TC, validate=TRUE))

expect_matrix_equal <- function(x, y, tol=1e-10, label="A", expected.label="B") {
  expect_equal(dim(x), dim(y), label=label, expected.label=expected.label)
  x <- as.matrix(x) # expect_equivalent doesn't compare sparse matrices
  y <- as.matrix(y) # and expect_equal doesn't seem to work because of method dispatch for compare()
  expect_equivalent(x, y, tolerance=tol, label=label, expected.label=expected.label)
}

## helper function: compare all relevant information from two DSM objects
compare.dsm <- function (A, B, check.marginals=TRUE, label="<TODO>") {
  name.A <- deparse(substitute(A))
  expect_failure(expect_error(check.dsm(A, validate=TRUE)))
  expect_failure(expect_error(check.dsm(B, validate=TRUE)))
  expect_matrix_equal(A$M, B$M, label=sprintf("%s$M", name.A), expected.label=sprintf("reference (%s)", label))
  expect_equal(A$rows$term, B$rows$term, label=sprintf("%s$rows$term", name.A), expected.label=sprintf("reference (%s)", label))
  expect_equal(A$cols$term, B$cols$term, label=sprintf("%s$cols$term", name.A), expected.label=sprintf("reference (%s)", label))
  expect_equal(A$rows$nnzero, B$rows$nnzero, label=sprintf("%s$rows$nnzero", name.A), expected.label=sprintf("reference (%s)", label))
  expect_equal(A$cols$nnzero, B$cols$nnzero, label=sprintf("%s$cols$nnzero", name.A), expected.label=sprintf("reference (%s)", label))
  if (check.marginals) {
    expect_equal(A$rows$f, B$rows$f, label=sprintf("%s$rows$f", name.A), expected.label=sprintf("reference (%s)", label))
    expect_equal(A$cols$f, B$cols$f, label=sprintf("%s$cols$f", name.A), expected.label=sprintf("reference (%s)", label))
    expect_equal(A$globals$N, B$globals$N, label=sprintf("%s$globals$N", name.A), expected.label=sprintf("reference (%s)", label))
  }
}


## Co-occurrence matrix in triplet representation
triplet.file <- system.file("extdata", "term_context_triplets.gz", package="wordspace", mustWork=TRUE)
test_that("triplet format can be read", {
  TC.triplet <- read.dsm.triplet(triplet.file, freq=TRUE, sort=TRUE, encoding="ascii")
  compare.dsm(TC.triplet, TC, label="DSM loaded from triplet file", check.marginals=FALSE) # correct marginal frequencies cannot be obtained from triplet file
})
  
## Triplet file with explicit marginal frequencies
marginal.file <- system.file("extdata", "term_context_marginals.txt.gz", package="wordspace", mustWork=TRUE)
test_that("triplet format with separate marginal frequencies can be read", {
  TC.triplet <- read.dsm.triplet(triplet.file, rowinfo=marginal.file, colinfo=marginal.file, N=TC$globals$N, freq=TRUE, sort=TRUE, encoding="ascii")
  compare.dsm(TC.triplet, TC, label="DSM loaded from triplet file with explicit marginals")
  expect_equal(TC.triplet$rows$df, TC.triplet$rows$nnzero)
  expect_equal(TC.triplet$cols$df, TC.triplet$cols$nnzero)
})

## span size adjustment
test_that("triplet format can be loaded with span size adjustment", {
  TC.spansize <- read.dsm.triplet(triplet.file, rowinfo=marginal.file, colinfo=marginal.file, N=TC$globals$N, freq=TRUE, span.size=42, sort=TRUE, encoding="ascii")
  expect_equal(TC.spansize$rows$f, TC$rows$f * 42)
  expect_equal(TC.spansize$cols$f, TC$cols$f)
})


## UCS export format
ucs.file <- system.file("extdata", "term_context_ucs", package="wordspace", mustWork=TRUE)
test_that("UCS export format can be read", {
  TC.ucs <- read.dsm.ucs(ucs.file, encoding="ascii")
  compare.dsm(TC.ucs, TC, label="DSM loaded from UCS export format")
})


## Check that non-ASCII characters in different encodings are read correctly
Encode.ref <- c("Test", "\u{0164}\u{00E9}\u{015F}t") # expected rownames of DSM

test.encoding <- function (filename, encoding, ref=Encode.ref) {
  test_that(sprintf("file encoding %s can be read correctly", encoding), {
    skip_if(!(encoding %in% iconvlist()), sprintf("character encoding %s not supported on this platform", encoding))
    model <- read.dsm.triplet(system.file("extdata", filename, package="wordspace", mustWork=TRUE), 
                              encoding=encoding, tokens=TRUE, sort=FALSE)
    expect_equal(model$rows$term, ref, label=filename, expected.label=sprintf("reference (%s)", encoding))
  })
}

test.encoding("tokens_utf8.txt", "UTF-8")
test.encoding("tokens_latin2.txt", "ISO-8859-2")
test.encoding("tokens_utf16.txt", "UTF-16")
