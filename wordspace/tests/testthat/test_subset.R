## Validate subset() methods for DSM objects
context("Subsetting DSMs")
library(wordspace)

## example from ?subset.dsm
model <- DSM_TermContext

test_that("subset() works for rows and columns", {
  m2 <- subset(model, nchar(term) <= 4, nnzero <= 3)
  expect_equal(dim(m2), c(3, 4))              # should extract a 3 x 4 matrix
  expect_equal(m2$rows$term, c("cat", "dog", "time")) # check that we have the expected rows
  expect_equal(m2$rows$term, rownames(m2$M))  # and that dimnames of the matrix are set correctly 
  expect_equal(m2$cols$nnzero, c(2, 2, 1, 0)) # check that nonzero counts have been updated
})

test_that("subset() drops zero rows and columns", {
  m2 <- subset(model, nnzero <= 3, nnzero <= 3, drop.zeroes=TRUE)
  expect_equal(dim(m2), c(3, 2)) # empty rows/columns should have been dropped
  expect_true(all(m2$rows$nnzero > 0))
  expect_true(all(m2$cols$nnzero > 0))
})

test_that("subset() can apply conditions recursively", {
  ## without recursive=TRUE, condition on nonzero count is only satisfied in original matrix
  m2 <- subset(model, nchar(term) <= 4, nnzero >= 2)
  expect_equal(dim(m2), c(3, 7))
  ## with recursive=TRUE, keep subsetting until both conditions are satisfied
  m3 <- subset(model, nchar(term) <= 4, nnzero >= 2, recursive=TRUE)
  expect_equal(dim(m3), c(3, 4))
})
