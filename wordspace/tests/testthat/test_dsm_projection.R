## Validate SVD and other latent subspace projections
## based on small "hieroglyphs" example matrix
context("Dimensionality reduction")
library(Matrix)
library(wordspace)

expect_matrix_equal <- function(x, y, tol=1e-10, ignore.sign=FALSE, note="", note1=note, note2=note) {
  name.x <- deparse(substitute(x))
  name.y <- deparse(substitute(y))
  x <- as.matrix(x)
  y <- as.matrix(y)
  expect_equal(dim(x), dim(y), 
               label=sprintf("dim(%s)%s", name.x, note1), 
               expected.label=sprintf("dim(%s)%s", name.y, note2))
  if (ignore.sign) {
    sign.vec <- sapply(seq_len(ncol(x)), function (i) {
      j <- which.max(abs(x[, i]))
      sign(x[j, i]) * sign(y[j, i])
    })
    x <- scaleMargins(x, cols=sign.vec)
  }
  expect_equivalent(x, y, tolerance=tol,
                    label=paste0(name.x, note1), 
                    expected.label=paste0(name.y, note2))
}

## compare distances between row vectors (for M1 and M2 in different subspaces)
##  - tol = proportion of average distance (for matrix with larger distances)
dist.compare <- function (M1, M2, method="euclidean", tol=.001, 
                          note="", note1=note, note2=note, expect=TRUE) { 
  name.M1 <- deparse(substitute(M1))
  name.M2 <- deparse(substitute(M2))
  if (expect) expect_true(nrow(M1) == nrow(M2),
                          label=sprintf("nrow(%s)%s == nrow(%s)%s", name.M1, note1, name.M2, note2))
  nr <- nrow(M1)
  M1 <- dist.matrix(M1, method=method)
  M2 <- dist.matrix(M2, method=method)
  n <- nr * (nr - 1)
  sum1 <- sum(M1) / n
  sum2 <- sum(M2) / n
  d <- sqrt(sum((M1 - M2)^2) / n) # root mean squared difference
  rel.error <- d / max(sum1, sum2)
  if (expect) expect_lte(rel.error, !!tol,
                         label=sprintf("relative error for distances in %s%s vs. %s%s", name.M1, note1, name.M2, note2))
  invisible(rel.error)
}

## test data: sparse co-occurrence matrix in sparse and dense representation
O <- DSM_HieroglyphsMatrix
E <- outer(rowSums(O), colSums(O)) / sum(O)
M1 <- scale((O - E) / sqrt(E), center=TRUE, scale=FALSE)  # z-score (signed), centered (so SVD projection == PCA)
M2 <- as(M1, "dgCMatrix") # force sparse matrix representation to test sparse algorithms


## full-rank SVD decomposition
svd.R <- svd(M1, nu=6, nv=6)
proj.R <- svd.R$u %*% diag(svd.R$d) # projection into SVD space
approx.R <- proj.R %*% t(svd.R$v)   # matrix approximation

proj.ws <- dsm.projection(M1, method="svd", n=6, with.basis=TRUE)
approx.ws <- proj.ws %*% t(attr(proj.ws, "basis"))

test_that("full-rank SVD matches original matrix", {
  dist.compare(proj.R, M1)
  dist.compare(approx.R, M1)
  expect_matrix_equal(approx.R, M1, tol=1e-6)
  
  dist.compare(proj.ws, M1)
  dist.compare(approx.ws, M1)
  expect_matrix_equal(approx.ws, M2, tol=1e-6)
})


## low-rank SVD projection
svd.R <- svd(M1, nu=3, nv=3)
proj.R <- svd.R$u %*% diag(svd.R$d[1:3]) # projection into SVD space
proj.ws <- dsm.projection(M1, method="svd", n=3)
proj.ws.sparse <- dsm.projection(M2, method="svd", n=3)

test_that("low-rank SVD projection is reasonable approximation", {
  dist.compare(proj.R, proj.ws)
  dist.compare(proj.R, M1, tol=.2) # relative error for SVD approximation should be ca. 10% 
  dist.compare(proj.ws, M1, tol=.2)
  dist.compare(proj.ws.sparse, M2, tol=.2)
  dist.compare(proj.ws, proj.ws.sparse)
  
  expect_matrix_equal(proj.R, proj.ws, tol=1e-6, ignore.sign=TRUE) # must ignore arbitrary sign of SVD dims
  expect_matrix_equal(proj.R, proj.ws.sparse, tol=1e-6, ignore.sign=TRUE)
})


## verify power scaling for dense and sparse SVD
test_that("SVD power scaling works as expected", {
  proj.ws.p2 <- dsm.projection(M1, method="svd", n=3, power=2)        # dense
  expect_matrix_equal(proj.ws.p2, svd.R$u %*% diag(svd.R$d[1:3]^2), tol=1e-6, ignore.sign=TRUE)
  proj.ws.sparse.p2 <- dsm.projection(M2, method="svd", n=3, power=2) # sparse
  expect_matrix_equal(proj.ws.sparse.p2, proj.ws.p2, tol=1e-6, ignore.sign=TRUE)
  
  expect_matrix_equal(proj.ws.p2, scaleMargins(proj.ws, cols=attr(proj.ws, "sigma")), tol=1e-6, ignore.sign=TRUE) # post-hoc power scaling
  expect_matrix_equal(proj.ws.sparse.p2, scaleMargins(proj.ws.sparse, cols=attr(proj.ws.sparse, "sigma")), tol=1e-6, ignore.sign=TRUE)
  
  proj.ws.p0 <- dsm.projection(M1, method="svd", n=3, power=0)       # whitening
  expect_matrix_equal(proj.ws.p0, svd.R$u, tol=1e-6, ignore.sign=TRUE)
})  


## full-rank randomized SVD
test_that("full-rank rSVD is equivalent to SVD", {
  set.seed(42) # for predictable results
  rsvd.proj.ws <- dsm.projection(M1, method="rsvd", n=6, oversampling=2, with.basis=TRUE)
  rsvd.approx.ws <- rsvd.proj.ws %*% t(attr(rsvd.proj.ws, "basis"))
  set.seed(42)
  rsvd.proj.ws.sparse <- dsm.projection(M2, method="rsvd", n=6, oversampling=2, with.basis=TRUE)
  rsvd.approx.ws.sparse <- rsvd.proj.ws.sparse %*% t(attr(rsvd.proj.ws.sparse, "basis"))
  
  dist.compare(rsvd.proj.ws, M1)
  expect_matrix_equal(rsvd.approx.ws, M1, tol=1e-6)
  dist.compare(rsvd.proj.ws.sparse, M2)
  expect_matrix_equal(rsvd.approx.ws.sparse, M2, tol=1e-6)
})  


## randomized SVD projection (with intermediate full rank)
test_that("rSVD with full-rank oversampling is equivalent to SVD", {
  set.seed(42)
  rsvd.proj.ws <- dsm.projection(M1, method="rsvd", n=3, oversampling=3)
  set.seed(42)
  rsvd.proj.ws.sparse <- dsm.projection(M2, method="rsvd", n=3, oversampling=3)
  
  dist.compare(rsvd.proj.ws, M1, tol=.2)
  dist.compare(rsvd.proj.ws, proj.ws)
  dist.compare(rsvd.proj.ws.sparse, M2, tol=.2)
  dist.compare(rsvd.proj.ws.sparse, proj.ws.sparse)
})


## randomized SVD projection (with reduced intermediate rank)
test_that("rSVD with reduced oversampling approximates to SVD", {
  set.seed(42)
  rsvd.proj.ws <- dsm.projection(M1, method="rsvd", n=2, oversampling=2)
  set.seed(42)
  rsvd.proj.ws.sparse <- dsm.projection(M2, method="rsvd", n=2, oversampling=2)

  dist.compare(rsvd.proj.ws, M1, tol=.35) # relative error should be ca. 30%
  dist.compare(rsvd.proj.ws.sparse, M2, tol=.35)
  dist.compare(rsvd.proj.ws, rsvd.proj.ws.sparse) # should be identical
})


## random indexing (with completely filled random vectors)
relerr.dense <- sapply(1:100, function (k) {
  set.seed(k)
  ri.proj.ws <- dsm.projection(M1, method="ri", n=3, rate=1, verbose=FALSE)
  dist.compare(ri.proj.ws, M1, expect=FALSE)
})
relerr.sparse <- sapply(1:100, function (k) {
  set.seed(k)
  ri.proj.ws <- dsm.projection(M2, method="ri", n=3, rate=1, verbose=FALSE)
  dist.compare(ri.proj.ws, M2, expect=FALSE)
})

res <- t.test(relerr.dense, relerr.sparse)
if (FALSE) {
  cat(sprintf("\nRelative error for dense RI (n=3, wordspace):  %4.1f +/- %4.1f %%\n", 
              100 * mean(relerr.dense), 100 * sd(relerr.dense)))
  cat(sprintf("Relative error for sparse RI (n=3, wordspace): %4.1f +/- %4.1f %%\n", 
              100 * mean(relerr.sparse), 100 * sd(relerr.sparse)))
  cat(sprintf("T-test for difference in means:  p = %.3f\n", res$p.value))
}
  
test_that("Random Indexing produces reasonable approximations", {
  expect_lt(min(relerr.dense), .3)
  expect_lt(median(relerr.dense), .5)
  expect_lt(min(relerr.sparse), .3)
  expect_lt(median(relerr.sparse), .5)
  expect_gt(res$p.value, .01) # with completely filled random vectors, both methods should have similar performance 
})
