library(testthat)
library(rpf)

context("goodness of fit")

# traditional Pearson X^2 goodness of fit test
pearson.gof <- function(observed, expected, df) {
  x2 <- sum((observed - expected)^2/expected)
  if (missing(df)) {
    df <- (dim(observed)[1]-1) * (dim(observed)[2]-1)
  }
  pchisq(x2, df, lower.tail=FALSE)
}

ms <- function(observed, expected, draws) {
  draws * sum((observed - expected)^2)
}

mc.gof.test <- function(observed, expected) {
  crosstabTest(matrix(observed, nrow=1), matrix(c(expected), nrow=1), 10000)
}

geolen <- function(v) sqrt(sum(v * v))
angle <- function(v1, v2) acos(sum(v1 * v2) / (geolen(v1) * geolen(v2)))

test_that("mc", {
  set.seed(1)
  bins <- 8L
  v1 <- matrix(runif(bins), nrow=2)
  v2 <- matrix(runif(bins), nrow=2)
  v1 <- v1 / sum(v1)
  v2 <- v2 / sum(v2)
  draws <- 300  # ptw2011 accuracy is proportional to draws
  
  got <- data.frame(gap=seq(.5, .05, -.025))
  for (xx in 1:dim(got)[1]) {
    gap <- got[xx,'gap']
    v3 <- v1 * (1-gap) + v2 * gap
    got$angle[xx] <- angle(v1, v3)
    observed <- v1 * draws
    expected <- v3 * draws
    got$mc[xx] <- mc.gof.test(observed, expected)
    got$ptw[xx] <- ptw2011.gof.test(observed, expected)
    got$x2[xx] <- pearson.gof(observed, expected, bins-1)
  }

  expect_equal(got$mc, got$ptw, .01)
#  print(max(abs(got$x2 - got$mc)))
  expect_true(max(abs(got$x2 - got$mc)) > .15)
})

test_that("crazy1", {
  # can obtain less than 0 without check
  observed <- structure(c(0L, 0L, 0L, 0L, 0L, 158L, 35L, 51L, 22L, 40L), .Dim = c(5L,  2L))
  expected <- structure(c(78.8628, 5.0351, 2.835, 0.5783, 0.4605, 108.4319,  23.5408, 35.5775, 17.2235, 33.4546), .Dim = c(5L, 2L))
  expect_true(ptw2011.gof.test(observed, expected) >= 0)
})

drawRandomProportion <- function(expected) {
  total <- sum(expected)
  prob <- expected / total
  sim <- rep(NA, length(expected))
  rowSim <- sample.int(length(expected), size=total, prob=prob, replace=TRUE)
  sim <- tabulate(rowSim, length(expected))
  sim
}

drawRandomCrosstab <- function(expected) {
  rowTotals <- apply(expected, 1, sum)
  sim <- matrix(NA, nrow(expected), ncol(expected))
  for (rx in 1:nrow(expected)) {
    rowSim <- sample.int(ncol(expected), size=rowTotals[rx],
                         prob=c(expected[rx,] / rowTotals[rx]), replace=TRUE)
    sim[rx,] <- as.numeric(tabulate(rowSim, ncol(expected)))
  }
  sim
}

if (0) {
  E = structure(c(2.849, 13.056, 19.17, 28.498, 41.552, 41.064, 56.528,  65.796, 62.687, 66.283, 56.303, 38.329, 49.191, 40.029, 31.637,  19.376, 11.155, 5.057, 2.257, 0.347, 0.151, 0.944, 1.83, 3.502,  6.448, 7.936, 13.472, 19.204, 22.313, 28.717, 29.697, 24.671,  38.809, 38.971, 38.363, 29.624, 21.845, 12.943, 7.743, 1.653), .Dim = c(20L,  2L))
  trials <- 150
  got <- rep(NA, trials)
  for (rep in 1:trials) {
    got[rep] <- crosstabTest(drawRandomCrosstab(E),E, 1000)
  }
  
  require(ggplot2)
  qplot(c(0, 1), stat = "function", fun = Vectorize(function(x) sum(got < x)/length(got)), geom = "line") +
    geom_abline(slope=1, color="red")+ coord_fixed()
}
