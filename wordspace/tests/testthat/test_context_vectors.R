## Test computation of centroids as context vectors
context("Context vectors")
library(wordspace)

vec.compare <- function (x, y, name="vector comparison", normalize=TRUE, tol=1e-10, verbose=TRUE) {
  x <- as.vector(x)
  y <- as.vector(y)
  if (any(is.na(x))) {
    cat(sprintf("%s: missing values in first vector\n", name))
    invisible(FALSE)
  } else if (any(is.na(y))) {
    cat(sprintf("%s: missing values in second vector\n", name))
    invisible(FALSE)
  } else if (length(x) != length(y)) {
    if (verbose) cat(sprintf("%s: different vector lengths %d != %d\n", name, length(x), length(y)))
    invisible(FALSE)
  } else {
    if (normalize) {
      x <- x / sqrt(sum(x ^ 2))
      y <- y / sqrt(sum(y ^ 2))
    }
    max.diff <- max(abs(x - y))
    if (max.diff < tol) {
      invisible(TRUE)
    } else {
      if (verbose) cat(sprintf("%s: largest difference between vectors = %g exceeds tolerance limit\n", name, max.diff))
      invisible(FALSE)
    }
  }
}

## test case: centroid vector for the document "cat cat dog cause cause"
M <- DSM_TermTermMatrix
x.ref <- (2 * M["cat", ] + M["dog", ] + 2 * M["cause", ]) / 5 

test_that("context.vectors() works with different context representations", {
  doc1 <- "cat cat dog cause cause" # centroid of document as string
  x1 <- context.vectors(M, doc1)
  expect_equivalent(x1, x.ref, tolerance=1e-12)
  doc2 <- c("cat", "cat", "dog", "cause", "cause") # centroid of document as list of tokens
  x2 <- context.vectors(M, list(doc2))
  expect_equivalent(x2, x.ref, tolerance=1e-12)
  doc3 <- c("cat", "cause", "cause", "dog", "cat") # should be independent of ordering
  x3 <- context.vectors(M, list(doc3))
  expect_equivalent(x3, x.ref, tolerance=1e-12)
  doc4 <- c(cat=1, cat=1, dog=1, cause=1, cause=1) # centroid of document as bag of words
  x4 <- context.vectors(M, list(doc4))
  expect_equivalent(x4, x.ref, tolerance=1e-12)
  doc5 <- c(cat=2, dog=1, cause=2) # with aggregated frequency counts
  x5 <- context.vectors(M, list(doc5))
  expect_equivalent(x5, x.ref, tolerance=1e-12)
  doc6 <- c(cat=.2, dog=.1, cause=.2) # check that non-integer values work correctly
  x6 <- context.vectors(M, list(doc6))
  expect_equivalent(x6, x.ref, tolerance=1e-12) # centroid should be independent of scaling of weights
  expect_equivalent(x5, x6, tolerance=1e-12) 
})
