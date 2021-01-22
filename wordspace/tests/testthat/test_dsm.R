## Test whether DSM object has expected structure
context("DSM object structure")
library(wordspace)

## Re-construct built-in DSM objects to verify they have the expected structure and attributes
reconstruct.dsm <- function (Orig) {
  ## extract relevant information in a "minimally similar" way
  if (dsm.is.canonical(Orig$M)$sparse) {
    M <- Orig$M
    attr(M, "nonneg") <- FALSE
  } else {
    M <- matrix(Orig$M, nrow(Orig$M), ncol(Orig$M), dimnames=dimnames(Orig$M))
  }
  row.f <- Orig$rows[sample(nrow(M)), c("term", "f")] # ordering of data frame shouldn't matter
  col.f <- Orig$cols[sample(ncol(M)), c("term", "f")]
  N <- Orig$globals$N
  dsm(M=M, rowinfo=row.f, colinfo=col.f, N=N, raw.freq=TRUE)
}

test_that("built-in DSM objects have expected structure", {
  expect_equal(DSM_TermTerm, reconstruct.dsm(DSM_TermTerm))
  expect_equal(DSM_TermContext, reconstruct.dsm(DSM_TermContext))
})
