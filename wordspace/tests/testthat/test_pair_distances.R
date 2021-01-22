## Validate pair.distances() implementation
context("Pair distances")
library(wordspace)

qw <- function (x) unlist(strsplit(x, "\\s+", perl=TRUE)) # Perl's qw()

M <- DSM_TermContextMatrix                # Manhattan distances are easy to work out by hand for this matrix
DM <- dist.matrix(M, method="manhattan")  # will be used to test direct lookup in distance matrix

## test data (including pairs with missing items)
w1 <-  qw("cat  cat  dog  rule  dog    rule")
w2 <-  qw("dog  pig  cat  cat   cause  pig")
d.exp <- c(24, Inf,  24,  Inf,  32,    Inf) # expected distances,
r.exp <- c(2,  Inf,  1,   Inf,  6,     Inf) # forward ranks,
b.exp <- c(1,  Inf,  2,   Inf,  4,     Inf) # and backward ranks

## compute pair distances from vectors
test_that("Manhattan distances for word pairs are correct", {
  d1 <- pair.distances(w1, w2, M, method="manhattan")
  expect_equivalent(d1, d.exp)
  expect_equivalent(names(d1), paste(w1, w2, sep="/")) # check correct labels
})

test_that("Manhattan ranks for word pairs are correct", {
  r1 <- pair.distances(w1, w2, M, method="manhattan", rank="fwd")
  expect_equivalent(r1, r.exp)
  b1 <- pair.distances(w1, w2, M, method="manhattan", rank="bwd")
  expect_equivalent(b1, b.exp)
})
  
## compute pair distances by direct lookup in pre-computed matrix
test_that("Manhattan distances can be looked up in pre-computed matrix", {
  d2 <- pair.distances(w1, w2, DM, method="euclidean") # verify that method will be ignored!
  expect_equivalent(d2, d.exp)
  expect_equivalent(names(d2), paste(w1, w2, sep="/"))
})
  
test_that("Manhattan ranks can be looked up in pre-computed matrix", {
  r2 <- pair.distances(w1, w2, DM, method="euclidean", rank="fwd")
  expect_equivalent(r2, r.exp)
  b2 <- pair.distances(w1, w2, DM, method="euclidean", rank="bwd")
  expect_equivalent(b2, b.exp)
})


## test data 2: looking up co-occurrence frequencies in a fake dist.matrix
DM <- as.distmat(DSM_TermContextMatrix, similarity=TRUE)
w1 <-  qw("cat    cat   time  time       yyz   cause      yyz")
w2 <-  qw("Feral  Rush  Boat  Back_Pain  Kant  Back_Pain  Rush")
d.exp <- c(7,     -Inf, 2,    0,         -Inf, 6,         -Inf)
r.exp <- c(3,     Inf,  1,    Inf,       Inf,  1,         Inf)
b.exp <- c(2,     Inf,  2,    Inf,       Inf,  1,         Inf)

## look up co-occurrence counts
test_that("pair.distances() can look up co-occurrence counts", {
  d3 <- pair.distances(w1, w2, DM)
  expect_equivalent(d3, d.exp)
  expect_equivalent(names(d3), paste(w1, w2, sep="/"))
  expect_true(attr(d3, "similarity")) # should be marked as similarity values
})

test_that("pair.distances() can compute neighbour ranks", {
  r3 <- pair.distances(w1, w2, DM, rank="fwd")
  expect_equivalent(r3, r.exp)
  b3 <- pair.distances(w1, w2, DM, rank="bwd")
  expect_equivalent(b3, b.exp)
})
