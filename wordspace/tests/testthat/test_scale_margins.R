## Validate scaleMargins() implementation against pure R code
## based on small "hieroglyphs" example matrix
context("Row and column scaling")
library(wordspace)
library(Matrix)

expect_matrix_equal <- function(x, y, tol=1e-6, note="", note1=note, note2=note) {
  name.x <- deparse(substitute(x))
  name.y <- deparse(substitute(y))
  expect_equal(dim(x), dim(y), 
               label=sprintf("dim(%s)%s", name.x, note1), 
               expected.label=sprintf("dim(%s)%s", name.y, note2))
  x <- as.matrix(x) # expect_equivalent doesn't compare sparse matrices
  y <- as.matrix(y) # (because all data items are stored in attributes)
  expect_equivalent(x, y, tolerance=tol,
                    label=paste0(name.x, note1), 
                    expected.label=paste0(name.y, note2))
}

## test data: sparse and dense co-occurrence matrix
M1 <- DSM_HieroglyphsMatrix  # dense matrix
M2 <- as(M1, "dgCMatrix")    # canonical sparse matrix

row.weights <- (1:nrow(M1)) - 1
col.weights <- (1:ncol(M1)) * 10

## pure R implementation as reference
scale.rows <- function (M, w) {
  if (length(w) == 1) w <- rep(w, nrow(M))
  stopifnot(length(w) == nrow(M))
  sweep(M, 1, w, FUN="*")
}
scale.cols <- function (M, w) {
  if (length(w) == 1) w <- rep(w, ncol(M))
  stopifnot(length(w) == ncol(M))
  sweep(M, 2, w, FUN="*")
}


## scale rows only (weight vector or constant)
test_that("rows of dense (M1) and sparse (M2) matrix can be scaled with individual weights", {
  M1.row.R <- scale.rows(M1, row.weights)
  M2.row.R <- scale.rows(M2, row.weights)
  M1.row.ws <- scaleMargins(M1, rows=row.weights)
  M2.row.ws <- scaleMargins(M2, rows=row.weights)
  
  expect_matrix_equal(M1.row.R, M2.row.R)
  expect_matrix_equal(M1.row.ws, M1.row.R)
  expect_matrix_equal(M2.row.ws, M2.row.R)
  expect_false(dsm.is.canonical(M1.row.ws)$sparse)
  expect_true(dsm.is.canonical(M2.row.ws)$sparse)
})

test_that("rows of dense (M1) and sparse (M2) matrix can be scaled with constant", {
  M1.row.R <- scale.rows(M1, 42)
  M2.row.R <- scale.rows(M2, 42)
  M1.row.ws <- scaleMargins(M1, rows=42)
  M2.row.ws <- scaleMargins(M2, rows=42)
  
  expect_matrix_equal(M1.row.R, M2.row.R)
  expect_matrix_equal(M1.row.ws, M1.row.R)
  expect_matrix_equal(M2.row.ws, M2.row.R)
  expect_false(dsm.is.canonical(M1.row.ws)$sparse)
  expect_true(dsm.is.canonical(M2.row.ws)$sparse)
})


## scale columns only (weight vector or constant)
test_that("columns of dense (M1) and sparse (M2) matrix can be scaled with individual weights", {
  M1.col.R <- scale.cols(M1, col.weights)
  M2.col.R <- scale.cols(M2, col.weights)
  M1.col.ws <- scaleMargins(M1, cols=col.weights)
  M2.col.ws <- scaleMargins(M2, cols=col.weights)

  expect_matrix_equal(M1.col.R,  M2.col.R)
  expect_matrix_equal(M1.col.ws, M1.col.R)
  expect_matrix_equal(M2.col.ws, M2.col.R)
  expect_false(dsm.is.canonical(M1.col.ws)$sparse)
  expect_true(dsm.is.canonical(M2.col.ws)$sparse)
})

test_that("columns of dense (M1) and sparse (M2) matrix can be scaled with constant", {
  M1.col.R <- scale.cols(M1, 42)
  M2.col.R <- scale.cols(M2, 42)
  M1.col.ws <- scaleMargins(M1, cols=42)
  M2.col.ws <- scaleMargins(M2, cols=42)
  
  expect_matrix_equal(M1.col.R,  M2.col.R)
  expect_matrix_equal(M1.col.ws, M1.col.R)
  expect_matrix_equal(M2.col.ws, M2.col.R)
  expect_false(dsm.is.canonical(M1.col.ws)$sparse)
  expect_true(dsm.is.canonical(M2.col.ws)$sparse)
})


## scale both rows and columns (weight vector or constant)
test_that("joint scaling of rows and columns works with individual weights", {
  M1.both.R <- scale.rows(scale.cols(M1, col.weights), row.weights)
  M2.both.R <- scale.rows(scale.cols(M2, col.weights), row.weights)
  M1.both.ws <- scaleMargins(M1, rows=row.weights, cols=col.weights)
  M2.both.ws <- scaleMargins(M2, rows=row.weights, cols=col.weights)

  expect_matrix_equal(M1.both.R,  M2.both.R)
  expect_matrix_equal(M1.both.ws, M1.both.R)
  expect_matrix_equal(M2.both.ws, M2.both.R)
  expect_false(dsm.is.canonical(M1.both.ws)$sparse)
  expect_true(dsm.is.canonical(M2.both.ws)$sparse)
})

test_that("joint scaling of rows and columns works with individual weights / constant", {
  M1.both.R <- scale.rows(scale.cols(M1, 42), row.weights)
  M2.both.R <- scale.rows(scale.cols(M2, 42), row.weights)
  M1.both.ws <- scaleMargins(M1, rows=row.weights, cols=42)
  M2.both.ws <- scaleMargins(M2, rows=row.weights, cols=42)

  expect_matrix_equal(M1.both.R,  M2.both.R)
  expect_matrix_equal(M1.both.ws, M1.both.R)
  expect_matrix_equal(M2.both.ws, M2.both.R)
  expect_false(dsm.is.canonical(M1.both.ws)$sparse)
  expect_true(dsm.is.canonical(M2.both.ws)$sparse)
})

test_that("joint scaling of rows and columns works with constant / individual weights", {
  M1.both.R <- scale.rows(scale.cols(M1, col.weights), 42)
  M2.both.R <- scale.rows(scale.cols(M2, col.weights), 42)
  M1.both.ws <- scaleMargins(M1, rows=42, cols=col.weights)
  M2.both.ws <- scaleMargins(M2, rows=42, cols=col.weights)

  expect_matrix_equal(M1.both.R,  M2.both.R)
  expect_matrix_equal(M1.both.ws, M1.both.R)
  expect_matrix_equal(M2.both.ws, M2.both.R)
  expect_false(dsm.is.canonical(M1.both.ws)$sparse)
  expect_true(dsm.is.canonical(M2.both.ws)$sparse)
})

test_that("joint scaling of rows and columns works with constants", {
  M1.both.R <- scale.rows(scale.cols(M1, 1/666), 42)
  M2.both.R <- scale.rows(scale.cols(M2, 1/666), 42)
  M1.both.ws <- scaleMargins(M1, rows=42, cols=1/666)
  M2.both.ws <- scaleMargins(M2, rows=42, cols=1/666)

  expect_matrix_equal(M1.both.R,  M2.both.R)
  expect_matrix_equal(M1.both.ws, M1.both.R)
  expect_matrix_equal(M2.both.ws, M2.both.R)
  expect_false(dsm.is.canonical(M1.both.ws)$sparse)
  expect_true(dsm.is.canonical(M2.both.ws)$sparse)
})


## test in-place operation (for internal use only)
test_that("rows and columns can be scaled with in-place operations", {
  M1.orig <- M1 + 0.0 # this should make copies
  M2.orig <- M2 + 0.0

  M1.both.R <- scale.rows(scale.cols(M1, col.weights), row.weights) # normal operation (returns new matrix)
  M2.both.R <- scale.rows(scale.cols(M2, col.weights), row.weights)
  M1.both.ws <- scaleMargins(M1, rows=row.weights, cols=col.weights)
  M2.both.ws <- scaleMargins(M2, rows=row.weights, cols=col.weights)

  expect_matrix_equal(M1.both.ws, M1.both.R)
  expect_matrix_equal(M2.both.ws, M2.both.R)
  expect_matrix_equal(M1, M1.orig) # with default duplicate=TRUE, copies of M1 and M2 are returned
  expect_matrix_equal(M2, M2.orig)

  M1.both.ws <- scaleMargins(M1, rows=row.weights, cols=col.weights, duplicate=FALSE) # in-place operation
  M2.both.ws <- scaleMargins(M2, rows=row.weights, cols=col.weights, duplicate=FALSE)

  expect_matrix_equal(M1.both.ws, M1.both.R)
  expect_matrix_equal(M2.both.ws, M2.both.R)
  expect_matrix_equal(M1, M1.both.R) # now M1 and M2 have been modified in-place
  expect_matrix_equal(M2, M2.both.R)
  ## in fact, DSM_HieroglpyhsMatrix will also be affected (since M1 is a reference to the same object)
  ## so make sure this is always the last test involving this matrix
})


## a sightly inefficient way to compute outer products
test_that("scaleMargins() corresponds to outer product", {
  x <- 1:5
  y <- c(0.25, 1, 2)
  outer.R <- outer(y, x)
  
  outer.ws <- matrix(1, nrow=length(y), ncol=length(x))
  outer.ws <- scaleMargins(outer.ws, cols=x, rows=y, duplicate=FALSE) # avoid extra copy
  
  expect_matrix_equal(outer.ws, outer.R)
})
