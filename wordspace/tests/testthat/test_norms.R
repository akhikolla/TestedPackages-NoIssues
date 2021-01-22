## Validate various matrix norms and normalization functions.
context("Vector norms")
library(wordspace)
library(Matrix)

## test data: small sparse and dense example matrix (unscaled)
M1 <- DSM_TermTermMatrix    # dense
M2 <- DSM_TermContextMatrix # sparse


## Euclidean, Manhattan and maximum norm, as well as Hamming length
test_that("Euclidean norm is correct for dense matrix", {
  tmp <- M1^2
  expect_equal(rowNorms(M1, "euclidean"), sqrt(rowSums(tmp)))
  expect_equal(colNorms(M1, "euclidean"), sqrt(colSums(tmp)))
})

test_that("Manhattan and maximum norms are correct for dense matrix", {
  tmp <- abs(M1)
  expect_equal(rowNorms(M1, "manhattan"), rowSums(tmp))
  expect_equal(colNorms(M1, "manhattan"), colSums(tmp))
  expect_equal(rowNorms(M1, "maximum"), apply(tmp, 1, max))
  expect_equal(colNorms(M1, "maximum"), apply(tmp, 2, max))
  expect_equal(rowNorms(M1, "minkowski", p=Inf), apply(tmp, 1, max))
  expect_equal(colNorms(M1, "minkowski", p=Inf), apply(tmp, 2, max))
})

test_that("Hamming length is correct for dense matrix", {
  tmp <- M1 != 0
  expect_equal(rowNorms(M1, "minkowski", p=0), rowSums(tmp))
  expect_equal(colNorms(M1, "minkowski", p=0), colSums(tmp))
})

test_that("Euclidean norm is correct for sparse matrix", {
  tmp <- M2^2
  expect_equal(rowNorms(M2, "euclidean"), sqrt(rowSums(tmp)))
  expect_equal(colNorms(M2, "euclidean"), sqrt(colSums(tmp)))
})

test_that("Manhattan and maximum norms are correct for sparse matrix", {
  tmp <- abs(M2)
  expect_equal(rowNorms(M2, "manhattan"), rowSums(tmp))
  expect_equal(colNorms(M2, "manhattan"), colSums(tmp))
  expect_equal(rowNorms(M2, "maximum"), apply(tmp, 1, max))
  expect_equal(colNorms(M2, "maximum"), apply(tmp, 2, max))
  expect_equal(rowNorms(M2, "minkowski", p=Inf), apply(tmp, 1, max))
  expect_equal(colNorms(M2, "minkowski", p=Inf), apply(tmp, 2, max))
})

test_that("Hamming length is correct for sparse matrix", {
  tmp <- M2 != 0
  expect_equal(rowNorms(M2, "minkowski", p=0), rowSums(tmp))
  expect_equal(colNorms(M2, "minkowski", p=0), colSums(tmp))
})


## various Minkowski norms
test_that("Minkowski norms are correct for dense and sparse matrix", {
  for (p in c(.1, .2, .5, 1, 1.5, 2, 3, 5, 10)) {
    q <- min(1 / p, 1)
    
    tmp <- abs(M1) ^ p
    expect_equal(rowNorms(M1, "minkowski", p=!!p), rowSums(tmp) ^ q)
    expect_equal(colNorms(M1, "minkowski", p=!!p), colSums(tmp) ^ q)
  
    tmp <- abs(M2) ^ p
    expect_equal(rowNorms(M2, "minkowski", p=!!p), rowSums(tmp) ^ q)
    expect_equal(colNorms(M2, "minkowski", p=!!p), colSums(tmp) ^ q)
  }
})
  

## validate row/column normalization (norms must be == 1)
test_that("rows and columns can be normalized", {
  nR1 <- nrow(M1); nC1 <- ncol(M1)
  nR2 <- nrow(M2); nC2 <- ncol(M2)
  
  for (norm in c("euclidean", "manhattan", "maximum")) {
    expect_equivalent(rowNorms(normalize.rows(M1, method=!!norm), method=!!norm), rep(1, nR1))
    expect_equivalent(colNorms(normalize.cols(M1, method=!!norm), method=!!norm), rep(1, nC1))
    expect_equivalent(rowNorms(normalize.rows(M2, method=!!norm), method=!!norm), rep(1, nR2))
    expect_equivalent(colNorms(normalize.cols(M2, method=!!norm), method=!!norm), rep(1, nC2))
  }

  for (p in c(.1, .2, .5, 1, 1.5, 2, 3, 5, 10)) {
    expect_equivalent(colNorms(normalize.cols(M1, method="minkowski", p=!!p), method="minkowski", p=!!p), rep(1, nR1)) 
    expect_equivalent(rowNorms(normalize.rows(M1, method="minkowski", p=!!p), method="minkowski", p=!!p), rep(1, nC1))
    expect_equivalent(rowNorms(normalize.rows(M2, method="minkowski", p=!!p), method="minkowski", p=!!p), rep(1, nR2))
    expect_equivalent(colNorms(normalize.cols(M2, method="minkowski", p=!!p), method="minkowski", p=!!p), rep(1, nC2))
  }
})

## error conditions: normalization not possible / reliable
test_that("normalization is not possible for Minkowski length with small p", {
  expect_error(normalize.rows(M1, method="minkowski", p=0))
  expect_error(normalize.cols(M1, method="minkowski", p=0))
  expect_error(normalize.rows(M1, method="minkowski", p=.01))
  expect_error(normalize.cols(M1, method="minkowski", p=.01))
})

## near-zero rows/columns are not normalized, but set to zero instead
fac <- c(1, 0, 1, 1, 1e-9, 0, 1e-12)
M1a <- scaleMargins(M1, rows=fac)
M1b <- scaleMargins(M1, cols=fac)
M2a <- scaleMargins(M2, rows=fac)
M2b <- scaleMargins(M2, cols=fac)

test_that("near-zero rows and columns are not normalized",{
  for (p in c(.5, .7, 1, 2, 5)) {
    tol <- if (p < 1) 0.1 else 1e-6 # scaling behaviour is very different for p < 1
    gold <- as.double(fac > 1e-6)
    M1a.norm <- normalize.rows(M1a, method="minkowski", p=p, tol=tol)
    M1b.norm <- normalize.cols(M1b, method="minkowski", p=p, tol=tol)
    M2a.norm <- normalize.rows(M2a, method="minkowski", p=p, tol=tol)
    M2b.norm <- normalize.cols(M2b, method="minkowski", p=p, tol=tol)
   
    
    expect_equivalent(rowNorms(M1a.norm, method="minkowski", p=!!p), gold)
    expect_equivalent(colNorms(M1b.norm, method="minkowski", p=!!p), gold)
    expect_equivalent(rowNorms(M2a.norm, method="minkowski", p=!!p), gold)
    expect_equivalent(colNorms(M2b.norm, method="minkowski", p=!!p), gold)
  }
})
