library(testthat)
library(rpf)

context("dataframe")

test_that("df basics", {
	df <- as.data.frame(matrix(sample.int(2, 5 * 100, replace=TRUE), 100, 5))
	cdf <- compressDataFrame(df)
	expect_true(nrow(cdf) < nrow(df))
	expect_equal(sum(cdf$freq), nrow(df))
	df2 <- expandDataFrame(cdf, "freq")
	expect_equal(nrow(df2), 100)
	expect_true(all(df2 == df[orderCompletely(df),]))
})

df <- as.data.frame(matrix(as.numeric(sample(2, 5 * 100, replace=TRUE)), 100, 5))
mask <- matrix(runif(5*100) < .1, ncol=5)
df[mask] <- NA
cdf <- compressDataFrame(df)
expect_equal(sum(cdf$freq), nrow(df))
df2 <- expandDataFrame(cdf, "freq")
expect_equal(nrow(df2), 100)
expect_true(all(df2 == df[orderCompletely(df),], na.rm=TRUE))
expect_true(all(is.na(df2) == is.na(df[orderCompletely(df),])))
