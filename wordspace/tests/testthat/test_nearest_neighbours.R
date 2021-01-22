## Validate nearest.neighbours() implementation
context("Nearest neighbours")
library(wordspace)

## example data from tutorial: hieroglyphs matrix
M <- log(DSM_HieroglyphsMatrix + 1) # similarity measure in tutorial is cosine on log frequencies

## nearest neighbours of mystery word "dog" as shown in tutorial
nn <- nearest.neighbours(M, "dog", convert=FALSE)
test_that("nearest neighbours of 'dog' correspond to tutorial", {
  expect_equal(names(nn), c("cat", "pig", "cup", "boat", "banana", "knife"))
  expect_equivalent(round(nn["cat"], 3), 0.961) # verify similarity values in tutorial slides
  expect_equivalent(round(nn["pig"], 3), 0.939)
  expect_equivalent(round(nn["knife"], 3), 0.770)
})

## extract same nearest neighbours with vector as target
test_that("nearest neighbours are found for vector", {
  nn2 <- nearest.neighbours(M, M["dog", ], convert=FALSE)
  expect_equal(nn2[-1], nn) # nn2 additionally contains "dog" as first NN
})

## return list of NN vectors for multiple targets
test_that("NN are found for multiple targets", {
  res <- nearest.neighbours(M, c("dog", "rat", "cat"), skip.missing=TRUE)
  expect_length(res, 2)
  expect_equal(names(res), c("dog", "cat"))
})

## return full distance matrix
DM1 <- nearest.neighbours(M, "dog", n=3, dist.matrix=TRUE)
attr(DM1, "selected") <- NULL # added by nearest.neighbours() function
DM2 <- dist.matrix(M, terms=c("dog", "cat", "pig", "cup"))
test_that("nearest.neighbours() can return distance matrix", {
  expect_equal(DM1, DM2) # nearest.neighbours() calls dist.matrix(), so results should be exactly the same
})
  
## access matrix by columns
test_that("NN can be computed for column vectors", {
  cM <- t(M)
  nn <- nearest.neighbours(cM, "dog", byrow=FALSE, convert=FALSE)
  expect_equal(names(nn), c("cat", "pig", "cup", "boat", "banana", "knife"))
  expect_equivalent(round(nn["cat"], 3), 0.961)
  cDM1 <- nearest.neighbours(cM, "dog", n=3, byrow=FALSE, dist.matrix=TRUE)
  attr(cDM1, "selected") <- NULL
  expect_equal(cDM1, DM2)
})


## find nearest neighbours in pre-computed distance matrix
pDM <- dist.matrix(M, method="cosine")
test_that("NN can be found in pre-computed distance matrix", {
  nn3 <- nearest.neighbours(pDM, "dog")
  expect_equal(names(nn3), names(nn))
  expect_equal(cos(nn3 * pi / 180), nn) # angular distance vs. cosine similarity
  nn4 <- nearest.neighbours(pDM, "dog", byrow=FALSE) # pDM is symmetric, so nn4 == nn3
  expect_equal(names(nn4), names(nn))
})

## search arbitrary similarity matrix for "nearest neighbours"
SMd <- as.distmat(DSM_TermTerm, similarity=TRUE)          # dense similarity matrix, extracted from DSM
SMs <- as.distmat(DSM_TermContextMatrix, similarity=TRUE) # sparse similarity matrix

test_that("nearest.neighbours() can be applied to arbitrary similarity matrix", {
  res1 <- nearest.neighbours(SMd, "time", n=3)
  expect_equal(res1, c(kill=134, likely=100, important=94))
  res2 <- nearest.neighbours(SMd, "feed", n=3, byrow=FALSE)
  expect_equal(res2, c(animal=86, dog=32, time=29))

  res3 <- nearest.neighbours(SMs, "animal", n=5) # should only return 3 nonempty cells
  expect_equal(res3, c(Pet=15, Feral=10, Felidae=2))
  res4 <- nearest.neighbours(SMs, "Back_Pain", n=3, byrow=FALSE) # only 2 nonempty cells
  expect_equal(res4, c(cause=6, reason=1))
})
