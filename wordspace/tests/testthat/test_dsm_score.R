## Validate C++ code in dsm_score() against pure R implementation
## based on small "hieroglyphs" example matrix
context("Scoring a DSM object")
library(wordspace)
library(Matrix)

## accept relatively large differences because C++ code might use different floating-point
## arithmetic than R code (depending on platform and optimizations)
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

## test data: co-occurrence matrix and corresponding DSM object
M <- DSM_HieroglyphsMatrix
M.dsm <- dsm(M, raw.freq=TRUE)


## marginals and expected frequencies
f1 <- rowSums(M)
f2 <- colSums(M)
N <- sum(M)
M.exp <- outer(f1, f2) / N

test_that("marginal and expected frequencies are correct", {
  expect_equivalent(f1, M.dsm$rows$f)
  expect_equivalent(f2, M.dsm$cols$f)
  expect_matrix_equal(N, M.dsm$globals$N)
})


## validate simple score functions and transformations
f <- M # raw co-occurrence frequency
f.dsm <- dsm.score(M.dsm, score="frequency", matrix.only=TRUE)
MI <- pmax(log2(M / M.exp), 0) # sparse pointwise Mutual Information
MI.dsm <- dsm.score(M.dsm, score="MI", matrix.only=TRUE)
t <- ifelse(M >= M.exp, (M - M.exp) / sqrt(M), 0) # sparse t-score
t.dsm <- dsm.score(M.dsm, score="t-score", matrix.only=TRUE)
tsqrt.dsm <- dsm.score(M.dsm, score="t-score", transform="root", matrix.only=TRUE) # square root transformation
tlog.dsm <- dsm.score(M.dsm, score="t-score", transform="log", matrix.only=TRUE) # logarithmic transformation

test_that("simple AMs and transformations are correct", {
  expect_matrix_equal(f, f.dsm)
  expect_matrix_equal(MI, MI.dsm)
  expect_matrix_equal(t, t.dsm)
  expect_matrix_equal(sqrt(t), tsqrt.dsm)
  expect_matrix_equal(log(t+1), tlog.dsm)
})


## validate score functions based on full contingency table
M.R1 <- outer(f1, rep(1, ncol(M)))
M.R2 <- N - M.R1
M.C1 <- outer(rep(1, nrow(M)), f2)
M.C2 <- N - M.C1
M.O11 <- M
M.O12 <- M.R1 - M.O11
M.O21 <- M.C1 - M.O11
M.O22 <- M.C2 - M.O12
M.E11 <- M.R1 * M.C1 / N # == M.exp
M.E12 <- M.R1 * M.C2 / N
M.E21 <- M.R2 * M.C1 / N
M.E22 <- M.R2 * M.C2 / N

G2 <- 2 * ( # log-likelihood (sparse and signed non-sparse
  ifelse(M.O11 > 0, M.O11 * log(M.O11 / M.E11), 0) +
  ifelse(M.O12 > 0, M.O12 * log(M.O12 / M.E12), 0) +
  ifelse(M.O21 > 0, M.O21 * log(M.O21 / M.E21), 0) +
  ifelse(M.O22 > 0, M.O22 * log(M.O22 / M.E22), 0)
)
G2.sparse <- ifelse(M.O11 > M.E11, G2, 0)
G2 <- ifelse(M.O11 >= M.E11, G2, -G2)
G2.sparse.dsm <- dsm.score(M.dsm, score="log-likelihood", matrix.only=TRUE)
G2.dsm <- dsm.score(M.dsm, score="log-likelihood", sparse=FALSE, negative.ok=TRUE, matrix.only=TRUE)

X2 <- # chi-squared (with Yates correction)
   N * (abs(M.O11 * M.O22 - M.O12 * M.O21) - N / 2) ^ 2 / (M.R1 * M.R2 * M.C1 * M.C2)
X2.sparse <- ifelse(M.O11 > M.E11, X2, 0)
X2 <- ifelse(M.O11 >= M.E11, X2, -X2)
X2.sparse.dsm <- dsm.score(M.dsm, score="chi-squared", matrix.only=TRUE)
X2.dsm <- dsm.score(M.dsm, score="chi-squared", sparse=FALSE, negative.ok=TRUE, matrix.only=TRUE)

test_that("complex AMs are correct", {
  expect_matrix_equal(G2.sparse, G2.sparse.dsm)
  expect_matrix_equal(G2, G2.dsm)
  expect_matrix_equal(X2.sparse, X2.sparse.dsm)
  expect_matrix_equal(X2, X2.dsm)
})


## verify that same results are obtianed on sparse matrix format
S <- as(M, "sparseMatrix")
stopifnot(dsm.is.canonical(S)$canonical && dsm.is.canonical(S)$sparse) # internal validation
S.dsm <- dsm(S, raw.freq=TRUE)

f.sparse <- dsm.score(S.dsm, score="frequency", matrix.only=TRUE)
MI.sparse <- dsm.score(S.dsm, score="MI", matrix.only=TRUE)
t.sparse <- dsm.score(S.dsm, score="t-score", matrix.only=TRUE)
tsqrt.sparse <- dsm.score(S.dsm, score="t-score", transform="root", matrix.only=TRUE)
tlog.sparse <- dsm.score(S.dsm, score="t-score", transform="log", matrix.only=TRUE)
G2.sparse.Sdsm <- dsm.score(S.dsm, score="log-likelihood", matrix.only=TRUE)
X2.sparse.Sdsm <- dsm.score(S.dsm, score="chi-squared", matrix.only=TRUE)

test_that("sparse AMs are identical on sparse matrix", {
  expect_is(f.sparse, "sparseMatrix")
  expect_matrix_equal(as.matrix(f.sparse), f.dsm)
  expect_is(MI.sparse, "sparseMatrix")
  expect_matrix_equal(as.matrix(MI.sparse), MI.dsm)
  expect_is(t.sparse, "sparseMatrix")
  expect_matrix_equal(as.matrix(t.sparse), t.dsm)
  expect_is(tsqrt.sparse, "sparseMatrix")
  expect_matrix_equal(as.matrix(tsqrt.sparse), tsqrt.dsm)
  expect_is(tlog.sparse, "sparseMatrix")
  expect_matrix_equal(as.matrix(tlog.sparse), tlog.dsm)
  expect_is(G2.sparse.Sdsm, "sparseMatrix")
  expect_matrix_equal(as.matrix(G2.sparse.Sdsm), G2.sparse)
  expect_is(X2.sparse.Sdsm, "sparseMatrix")
  expect_matrix_equal(as.matrix(X2.sparse.Sdsm), X2.sparse)
})


## pseudo-sparse AM (allowing negative scores, but only for nonzero entries of sparse matrix)
G2.nz.dsm <- dsm.score(M.dsm, score="log-likelihood", sparse=FALSE, negative.ok="nonzero", matrix.only=TRUE)
X2.nz.dsm <- dsm.score(M.dsm, score="chi-squared", sparse=FALSE, negative.ok="nonzero", matrix.only=TRUE)
G2.nz <- ifelse(M.O11 > 0, G2, 0) # for sparse matrix, compute signed scores only for nonzero cells
G2.nz.Sdsm <- dsm.score(S.dsm, score="log-likelihood", sparse=FALSE, negative.ok="nonzero", matrix.only=TRUE)
X2.nz <- ifelse(M.O11 > 0, X2, 0)
X2.nz.Sdsm <- dsm.score(S.dsm, score="chi-squared", sparse=FALSE, negative.ok="nonzero", matrix.only=TRUE)

test_that("pseudo-sparse AMs are computed correctly", {
  expect_matrix_equal(G2.nz.dsm, G2)
  expect_matrix_equal(X2.nz.dsm, X2)
  expect_matrix_equal(as.matrix(G2.nz.Sdsm), G2.nz)
  expect_matrix_equal(as.matrix(X2.nz.Sdsm), X2.nz)
})


## column scaling
tsqrt.std <- scale(sqrt(t)) # standardization
tsqrt.std.dsm <- dsm.score(M.dsm, score="t-score", transform="root", scale="standardize", matrix.only=TRUE) 

test_that("column standardization works", {
  expect_equivalent(apply(tsqrt.std.dsm, 2, mean), rep(0, ncol(M)))
  expect_equivalent(apply(tsqrt.std.dsm, 2, sd), rep(1, ncol(M)))
  expect_matrix_equal(tsqrt.std.dsm, tsqrt.std)
})

tsqrt.scale <- scale(sqrt(t), center=FALSE) # scale to unit variance without centering
tsqrt.scale.dsm <- dsm.score(M.dsm, score="t-score", transform="root", scale="scale", matrix.only=TRUE)
rms <- function (x) sqrt(sum(x * x) / (length(x) - 1)) # root mean square (cf. ?scale)

test_that("column scaling works", {
  expect_equivalent(apply(tsqrt.scale.dsm, 2, rms), rep(1, ncol(M)))
  expect_matrix_equal(tsqrt.scale.dsm, tsqrt.scale)
})
  
## NB: sparse DSM doesn't allow full standardization unless negative.ok=TRUE (and this is not recommended)
tsqrt.scale.sparse <- dsm.score(S.dsm, score="t-score", transform="root", scale="scale", matrix.only=TRUE) # scale to unit variance without centering

test_that("column scaling works for sparse matrix, too", {
  expect_equivalent(apply(tsqrt.scale.sparse, 2, rms), rep(1, ncol(M)))
  expect_matrix_equal(as.matrix(tsqrt.scale.sparse), tsqrt.scale)
})


## row normalization
tsqrt.norm <- sqrt(t) / rowSums(abs(sqrt(t))) # simplest case: Manhattan norm
tsqrt.norm.dsm <- dsm.score(M.dsm, score="t-score", transform="root", normalize=TRUE, method="manhattan", matrix.only=TRUE)
tsqrt.norm.sparse <- dsm.score(S.dsm, score="t-score", transform="root", normalize=TRUE, method="manhattan", matrix.only=TRUE)

test_that("row normalization works", {
  expect_matrix_equal(tsqrt.norm.dsm, tsqrt.norm)
  expect_matrix_equal(as.matrix(tsqrt.norm.sparse), tsqrt.norm)
})


## validate user-defined association measures against built-in C code
my.MI <- function (O11, E11, ...) log2(O11 / E11)
my.tscore <- function (O, E, ...) (O - E) / sqrt(O)
my.tscore.disc <- function (O, E, ...) (O - E) / sqrt(O + 1) # discounted t-score used as non-sparse measure
my.Dice <- function (O11, R1, C1, ...) 2 * O11 / (R1 + C1)
my.ll <- function (O, E, ...) {
  ll <- 2 * (ifelse(O > 0, O * log(O / E), 0) - (O - E))
  ifelse(O >= E, ll, -ll)
}
my.G2 <- function (O11, O12, O21, O22, E11, E12, E21, E22, ...) {
  ll <- ifelse(O11 > 0, O11 * log(O11 / E11), 0) +
    ifelse(O12 > 0, O12 * log(O12 / E12), 0) +
    ifelse(O21 > 0, O21 * log(O21 / E21), 0) +
    ifelse(O22 > 0, O22 * log(O22 / E22), 0)
  ifelse(O11 >= E11, 2 * ll, -2 * ll)
}
my.tfidf <- function (O11, cols, ..., K=7) O11 * log(1 / cols$df) # will need to annotate col$df explicitly

TT <- DSM_TermTerm    # dense co-occurrence matrix
TT$cols <- transform(TT$cols, df=((nnzero + 1) / (nrow(TT) + 1)))
TC <- DSM_TermContext # sparse co-occurrence matrix
TC$cols <- transform(TC$cols, df=((nnzero + 1) / (nrow(TC) + 1)))

test.AM <- function(model, AM.name, AM.fun, transform="none", sparse=TRUE, negative.ok=!sparse, force.dense=FALSE) {
  built.in <- dsm.score(model, AM.name, transform=transform, sparse=sparse, negative.ok=negative.ok, matrix.only=TRUE)
  user.defined <- dsm.score(model, AM.fun, transform=transform, sparse=sparse, negative.ok=negative.ok, matrix.only=TRUE, batchsize=10)
  if (force.dense) built.in <- as.matrix(built.in)
  label <- sprintf(" (%s%s%s)", 
                   if (transform == "none") "" else paste0(transform, " "),
                   AM.name,
                   if (sparse) if (negative.ok == "nonzero") ", pseudo-sparse" else ", sparse" else "")
  expect_equal(class(user.defined), class(built.in))
  expect_matrix_equal(user.defined, built.in, note2=label)
}

## dense matrix, sparse AM
test_that("user-defined AM match built-in AM (dense matrix, sparse AM)", {
  test.AM(TT, "MI", my.MI)
  test.AM(TT, "t-score", my.tscore)
  test.AM(TT, "Dice", my.Dice)
  test.AM(TT, "simple-ll", my.ll)
  test.AM(TT, "log-likelihood", my.G2)
  test.AM(TT, "tf.idf", my.tfidf)
  
  test.AM(TT, "simple-ll", my.ll, transform="log")
  test.AM(TT, "simple-ll", my.ll, transform="root")
  test.AM(TT, "simple-ll", my.ll, transform="sigmoid")
})

## dense matrix, dense AM
test_that("user-defined AM match built-in AM (dense matrix, dense AM)", {
  test.AM(TT, "MI", my.MI, sparse=FALSE)
  test.AM(TT, "t-score", my.tscore.disc, sparse=FALSE)
  test.AM(TT, "Dice", my.Dice, sparse=FALSE)
  test.AM(TT, "simple-ll", my.ll, sparse=FALSE)
  test.AM(TT, "log-likelihood", my.G2, sparse=FALSE)
  test.AM(TT, "tf.idf", my.tfidf, sparse=FALSE)

  test.AM(TT, "simple-ll", my.ll, transform="log", sparse=FALSE)
  test.AM(TT, "simple-ll", my.ll, transform="root", sparse=FALSE)
  test.AM(TT, "simple-ll", my.ll, transform="sigmoid", sparse=FALSE)
})

## sparse matrix, sparse AM
test_that("user-defined AM match built-in AM (sparse matrix, sparse AM)", {
  test.AM(TC, "MI", my.MI)
  test.AM(TC, "t-score", my.tscore)
  test.AM(TC, "Dice", my.Dice)
  test.AM(TC, "simple-ll", my.ll)
  test.AM(TT, "log-likelihood", my.G2)
  test.AM(TC, "tf.idf", my.tfidf)

  test.AM(TC, "simple-ll", my.ll, transform="log")
  test.AM(TC, "simple-ll", my.ll, transform="root")
  test.AM(TC, "simple-ll", my.ll, transform="sigmoid")
})

## sparse matrix, dense AM (casts to dense matrix)
test_that("user-defined AM match built-in AM (sparse matrix, dense AM)", {
  test.AM(TC, "MI", my.MI, sparse=FALSE)
  test.AM(TC, "t-score", my.tscore.disc, sparse=FALSE)
  test.AM(TC, "Dice", my.Dice, sparse=FALSE, force.dense=TRUE) # built-in AMs known to be sparse don't cast automatically
  test.AM(TC, "simple-ll", my.ll, sparse=FALSE)
  test.AM(TT, "log-likelihood", my.G2, sparse=FALSE)
  test.AM(TC, "tf.idf", my.tfidf, sparse=FALSE, force.dense=TRUE)
  
  test.AM(TC, "simple-ll", my.ll, transform="log", sparse=FALSE)
  test.AM(TC, "simple-ll", my.ll, transform="root", sparse=FALSE)
  test.AM(TC, "simple-ll", my.ll, transform="sigmoid", sparse=FALSE)
})

## sparse matrix, dense AM (only for nonzero cells)
test_that("user-defined AM match built-in AM (sparse matrix, quasi-sparse AM)", {
  test.AM(TC, "MI", my.MI, sparse=FALSE, negative.ok="nonzero")
  test.AM(TC, "t-score", my.tscore.disc, sparse=FALSE, negative.ok="nonzero")
  test.AM(TC, "Dice", my.Dice, sparse=FALSE, negative.ok="nonzero")
  test.AM(TC, "simple-ll", my.ll, sparse=FALSE, negative.ok="nonzero")
  test.AM(TT, "log-likelihood", my.G2, sparse=FALSE, negative.ok="nonzero")
  test.AM(TC, "tf.idf", my.tfidf, sparse=FALSE, negative.ok="nonzero")
  
  test.AM(TC, "simple-ll", my.ll, transform="log", sparse=FALSE, negative.ok="nonzero")
  test.AM(TC, "simple-ll", my.ll, transform="root", sparse=FALSE, negative.ok="nonzero")
  test.AM(TC, "simple-ll", my.ll, transform="sigmoid", sparse=FALSE, negative.ok="nonzero")
})
