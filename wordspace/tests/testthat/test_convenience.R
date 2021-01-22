## Test various convenience methods for DSM objects and related functionality
context("Convenience features")
library(wordspace)

dsm1 <- DSM_TermContext
dsm2 <- subset(DSM_TermTerm, f < 100000) # reduce to 5x7 matrix
dsm2 <- dsm.score(dsm2, function(O, E, ...) round(log2(O/E), 2))


## dim() of DSM object
test_that("dim() method works", {
  expect_equal(dim(dsm1), c(7, 7))
  expect_equal(nrow(dsm2), 5)
  expect_equal(ncol(dsm2), 7)
})


## reading and setting dimnames
rn <- c("cat", "dog", "animal", "reason", "cause") # expected row/column names
cn <- c("breed", "tail", "feed", "kill", "important", "explain", "likely")
test_that("dimnames() can be read correctly", {
  expect_equal(dimnames(dsm1), dimnames(dsm1$M))
  expect_equal(dimnames(dsm2), dimnames(dsm2$M))
  expect_equal(rownames(dsm2), rn)
  expect_equal(colnames(dsm2), cn)
})

test_that("dimnames() can be modified", {
  dsm2n <- dsm2 # modification of row/column names
  rownames(dsm2n)[3] <- "pet"
  expect_equal(dimnames(dsm2), list(rn, cn)) # must not affect original object
  expect_equal(colnames(dsm2n), cn) # no changes to column names yet
  rn2 <- rn; rn2[3] <- "pet" 
  expect_equal(rownames(dsm2n), rn2) # but one element of rownames has been changed
  cn2 <- LETTERS[1:7] # replace all column names
  colnames(dsm2n) <- cn2
  expect_equal(dimnames(dsm2), list(rn, cn))
  expect_equal(rownames(dsm2n), rn2)
  expect_equal(colnames(dsm2n), cn2)
  expect_failure(expect_error(check.dsm(dsm2n, validate=TRUE))) # checks that rows$term and cols$term have been updated so validation does not throw an error
})


## extraction of co-occurrence or score matrix
test_that("as.matrix() extracts appropriate co-occurrence matrix", {
  expect_equal(as.matrix(dsm1), dsm1$M)
  expect_equal(as.matrix(dsm2), dsm2$S) # automatic selection of S if available
  expect_equal(as.matrix(dsm2, what="M"), dsm2$M)
  expect_equal(as.matrix(dsm2, what="S"), dsm2$S)
})


## transposition
test_that("DSM object can be transposed", {
  dsm2t <- t(dsm2)
  expect_equal(nrow(dsm2t), ncol(dsm2))
  expect_equal(ncol(dsm2t), nrow(dsm2))
  expect_equal(dsm2t$M, t(dsm2$M))
  expect_equal(dsm2t$S, t(dsm2$S))
})


## efficient checks for non-negativity asnd nonzero count
## (validate against non-optimized pure R implementations)
R.signcount <- function (x) {
  if (is(x, "Matrix")) x <- as.matrix(x) # just test on small sparse matrix examples
  c(pos=sum(x > 0), zero=sum(x == 0), neg=sum(x < 0))
}
R.nonneg <- function (x) {
  if (is(x, "Matrix")) x <- as.matrix(x)
  !any(x < 0)
}
R.nnzero <- function (x) {
  if (is(x, "Matrix")) x <- as.matrix(x)
  sum(x != 0)
}

test_that("signcount() works for numeric and integer vectors", {
  x <- round(rnorm(1e6), .1) # test on a large numeric vector
  y <- as.integer(round(x))  # and integer counterpart
  expect_equal(signcount(x), R.signcount(x))
  expect_equal(sum(signcount(x)), length(x))
  expect_equal(signcount(y), R.signcount(y))
  expect_equal(sum(signcount(y)), length(y))
  expect_equal(signcount(x, "nonneg"), R.nonneg(x))
  expect_equal(signcount(y, "nonneg"), R.nonneg(y))
  expect_equal(signcount(x, "nnzero"), R.nnzero(x))
  expect_equal(signcount(y, "nnzero"), R.nnzero(y))
})

test_that("signcount() works for dense matrices", {
  expect_equal(signcount(DSM_TermTermMatrix), R.signcount(DSM_TermTermMatrix)) # dense numeric matrix
  M <- Matrix(DSM_HieroglyphsMatrix) # and a dense Matrix object
  expect_equal(signcount(M), R.signcount(M))
  expect_equal(sum(signcount(M)), prod(dim(M)))
  expect_equal(signcount(M, "nonneg"), R.nonneg(M))
  expect_equal(signcount(M, "nnzero"), R.nnzero(M))
})

test_that("signcount() works for sparse matrices", {
  expect_equal(signcount(DSM_TermContextMatrix), R.signcount(DSM_TermContextMatrix)) # sparse matrix (dgCMatrix)
  M <- as(DSM_TermContextMatrix, "dgTMatrix") # triplet representation is not supported
  expect_error(signcount(M), "must be.*vector")
  M <- as(as.matrix(M), "dgRMatrix") # sparse matrix (dgRMatrix), has only minimal support in Matrix package
  expect_equal(signcount(M), R.signcount(M))
  expect_equal(signcount(M, "nonneg"), R.nonneg(M))
  expect_equal(signcount(M, "nnzero"), R.nnzero(M))
})
