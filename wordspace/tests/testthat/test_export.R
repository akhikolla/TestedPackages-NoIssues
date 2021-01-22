## Test exporting and importing a DSM matrix in different file formats
context("Export/import functions")
library(wordspace)


## export word2vec text format
word2vec_hiero <- system.file("extdata", "word2vec_hiero.txt", package="wordspace", mustWork=TRUE)
M <- (DSM_HieroglyphsMatrix + .1) / 100 # make sure we don't depend on round half conventions
gold.lines <- readLines(word2vec_hiero)

test_that("word2vec format can be exported", {
  export_word2vec <- tempfile(fileext=".txt")
  write.dsm.matrix(M, export_word2vec, format="word2vec", round=1)
  expect_equal(readLines(export_word2vec), gold.lines)

  export_word2vec.gz <- tempfile(fileext=".txt.gz")
  fh <- gzfile(export_word2vec.gz, encoding="UTF-8")
  write.dsm.matrix(M, fh, format="word2vec", round=1)
  expect_equal(readLines(export_word2vec.gz), gold.lines)
})

## import word2vec text format
test_that("word2vec format can be imported", {
  M.word2vec <- read.dsm.matrix(word2vec_hiero, format="word2vec")
  colnames(M.word2vec) <- colnames(M) # no colnames in word2vec format
  expect_equal(M.word2vec, round(M, 1))
})

## round trip for larger matrix
BNC <- with(DSM_VerbNounTriples_BNC,
            dsm(target=noun, feature=verb, score=f, raw.freq=TRUE))
BNC <- dsm.score(BNC, score="simple-ll", transform="log", normalize=TRUE, update.nnzero=TRUE)
BNC100 <- dsm.projection(BNC, n=100, method="svd")

test_that("round trip to word2vec format works for larger matrix", {
  fn <- tempfile(fileext=".txt")
  write.dsm.matrix(BNC100, fn, format="word2vec", round=3)
  RoundTrip <- read.dsm.matrix(fn, format="word2vec")

  Gold <- round(BNC100[, ], 3) # round the original matrix,
  colnames(Gold) <- NULL       # removing extra attributes and column names
  expect_equal(RoundTrip, Gold)
})