library(testthat)
context("LOPART")

x <- c(1.1, 2.2, 5.5, 6.6)
no.labels <- data.frame(
  start=integer(),
  end=integer(),
  changes=integer())

test_that("LOPART with no labels and penalty=0", {
  out_list <- LOPART::LOPART(x, no.labels, 0)
  expect_equal(out_list$loss$changes_total, 3)
  expect_equal(out_list$loss$total_loss + sum(x*x), 0)
})

test_that("LOPART with no labels and big penalty", {
  penalty <- 100000
  out_list <- LOPART::LOPART(x, no.labels, penalty)
  expect_equal(out_list$loss$changes_total, 0)
  expect_equal(out_list$cost$mean, cumsum(x)/seq_along(x))
})
test_that("LOPART with no labels and penalty=Inf", {
  out_list <- LOPART::LOPART(x, no.labels, Inf)
  expect_equal(out_list$loss$changes_total, 0)
  expected.cost <- c(rep(Inf, 3), sum((mean(x)-x)^2-x^2))
  expect_equal(out_list$cost$cost_optimal, expected.cost)
})
test_that("LOPART_interface with no labels and penalty_unlabeled=Inf", {
  out_df <- LOPART::LOPART_interface(
    x, integer(), integer(), integer(),
    n_updates=length(x),
    penalty_unlabeled = Inf,
    penalty_labeled=Inf)
  expect_equal(out_df$mean, cumsum(x)/seq_along(x))
  n.changes <- sum(0 <= out_df$last_change)
  expect_equal(n.changes, 0)
})

test_that("LOPART with one positive label on [1,4] and penalty=0", {
  pos14.label <- data.frame(
    start=1,
    end=length(x),
    changes=1)
  out_list <- LOPART::LOPART(x, pos14.label, 0)
  expect_equal(out_list$loss$changes_total, 1)
  expect_equal(out_list$loss$penalized_cost, out_list$loss$total_loss)
})

test_that("LOPART with one positive label on [1,4] and penalty=Inf", {
  pos14.label <- data.frame(
    start=1,
    end=length(x),
    changes=1)
  out_list <- LOPART::LOPART(x, pos14.label, Inf)
  expect_equal(out_list$loss$changes_total, 1)
  expect_equal(out_list$loss$penalized_cost, Inf)
  m <- c(rep(mean(x[1:2]), 2), rep(mean(x[3:4]), 2))
  expected.loss <- sum((x-m)^2-x^2)
  expect_equal(out_list$loss$total_loss, expected.loss)
})

test_that("LOPART with one negative label on [1,4] and penalty=0", {
  neg14.label <- data.frame(
    start=1,
    end=length(x),
    changes=0)
  out_list <- LOPART::LOPART(x, neg14.label, 0)
  expect_equal(out_list$loss$changes_total, 0)
  m <- mean(x)
  expected.segs <- data.table::data.table(start=1L, end=4L, mean=m)
  expect_equal(out_list$segments, expected.segs)
  expected.cost <- sum((x-m)^2 - x^2)
  expect_equal(out_list$cost$cost_optimal, c(Inf, Inf, Inf, expected.cost))
})

pos13.label <- data.frame(
  start=1,
  end=3,
  changes=1)
test_that("LOPART with one positive label on [1,3] and penalty=0", {
  out_list <- LOPART::LOPART(x, pos13.label, 0)
  expect_equal(out_list$loss$changes_total, 2)
  expect_equal(out_list$segments$end, 2:4)
})
test_that("LOPART with one positive label on [1,3] and small penalty", {
  penalty <- 0.1
  out_list <- LOPART::LOPART(x, pos13.label, penalty)
  m <- c(rep(mean(x[1:2]), 2), x[3:4])
  expected.changes <- sum(diff(m) != 0)
  expected.loss <- sum((x-m)^2 - x^2)
  expect_equal(out_list$loss$changes_total, 2)
  expect_equal(out_list$loss$total_loss, expected.loss)
  expected.cost <- expected.loss+penalty*expected.changes
  expect_equal(out_list$loss$penalized_cost, expected.cost)
  expect_equal(out_list$segments$end, 2:4)
})
test_that("LOPART with one positive label on [1,3] and big penalty", {
  out_list <- LOPART::LOPART(x, pos13.label, 100000)
  expect_equal(out_list$segments$end, c(2, 4))
})
test_that("LOPART with one positive label on [1,3] and penalty=Inf", {
  out_list <- LOPART::LOPART(x, pos13.label, Inf)
  expect_equal(out_list$segments$end, c(2, 4))
})
test_that("LOPART_interface with one positive label on [1,3] and penalty=Inf", {
  out_df <- LOPART::LOPART_interface(
    x, pos13.label$start, pos13.label$end, pos13.label$changes,
    n_updates=length(x),
    penalty_labeled = 0,
    penalty_unlabeled = Inf)
  is.change <- 0 <= out_df$last_change
  change.vec <- out_df$last_change[is.change]
  end.vec <- c(change.vec+1, length(x))
  expect_equal(end.vec, c(2, 4))
})

two.labels <- data.frame(
  start=1:2,
  end=2:3,
  changes=1)
test_that("LOPART with two positive labels on [1,3] and big penalty", {
  out_list <- LOPART::LOPART(x, two.labels, 100000)
  expect_equal(out_list$segments$end, c(1, 2, 4))
})
test_that("LOPART with two positive labels on [1,3] and penalty=0", {
  out_list <- LOPART::LOPART(x, two.labels, 0)
  expect_equal(out_list$segments$end, 1:4)
})

three.labels <- data.frame(
  start=1:3,
  end=2:4,
  changes=c(1, 0, 1))
m <- c(x[1], rep(mean(x[2:3]), 2), x[4])
expected.loss <- sum((x-m)^2 - x^2)
test_that("LOPART with three labels and penalty=0", {
  out_list <- LOPART::LOPART(x, three.labels, 0)
  expect_equal(out_list$segments$end, c(1, 3, 4))
  expect_equal(out_list$loss$penalized_cost, expected.loss)
  expect_equal(out_list$loss$total_loss, expected.loss)
})
test_that("LOPART with three labels and big penalty", {
  penalty <- 100000
  out_list <- LOPART::LOPART(x, three.labels, penalty)
  expect_equal(out_list$segments$end, c(1, 3, 4))
  expect_equal(out_list$loss$penalized_cost, expected.loss+penalty*2)
  expect_equal(out_list$loss$total_loss, expected.loss)
})
test_that("error for negative penalty", {
  expect_error({
    LOPART::LOPART_interface(x, integer(), integer(), integer(), 1, -1)
  }, "penalty must be non-negative")
})

test_that("error for NA penalty", {
  expect_error({
    LOPART::LOPART_interface(x, integer(), integer(), integer(), 1, NA)
  }, "penalty must be non-negative")
})

test_that("error for label with start=end", {
  expect_error({
    LOPART::LOPART_interface(x, 1, 1, 1, 1, 1)
  }, "each label start must be less than its end")
})    

test_that("error for label with start>end", {
  expect_error({
    LOPART::LOPART_interface(x, 2, 1, 1, 1, 1)
  }, "each label start must be less than its end") 
})    

test_that("error for label changes not 0/1", {
  expect_error({
    LOPART::LOPART_interface(x, 1, 2, 5, 1, 1)
  }, "labeled number of changes must be 0 or 1")
})    

test_that("error for label start < prev label end", {
  expect_error({
    LOPART::LOPART_interface(x, c(1, 2), c(3, 4), c(1, 1), 1, 1)
  }, "each label start must be on or after previous end")
})

test_that("error for start/changes sizes that do not match", {
  expect_error({
    LOPART::LOPART_interface(x, c(1, 2), c(3, 4), 1, 1, 1)
  }, "input_label_start and input_label_changes sizes must match")
})

test_that("error for end/changes sizes that do not match", {
  expect_error({
    LOPART::LOPART_interface(x, 1, c(3, 4), 1, 1, 1)
  }, "input_label_end and input_label_changes sizes must match")
})

test_that("error for end equal to n data", {
  expect_error({
    LOPART::LOPART_interface(x, 1, 4, 1, 1, 1)
  }, "label end must be less than n data")
})

test_that("error for end > n data", {
  expect_error({
    LOPART::LOPART_interface(x, 1, 400, 1, 1, 1)
  }, "label end must be less than n data")
})

test_that("error for start < 0", {
  expect_error({
    LOPART::LOPART_interface(x, -1, 2, 1, 1, 1)
  }, "label start must be zero or larger")
})

test_that("error for no data", {
  expect_error({
    LOPART::LOPART_interface(numeric(), integer(), integer(), integer(), 1, 1)
  }, "no data")
})

test_that("error for NA data", {
  expect_error({
    LOPART::LOPART_interface(NA, integer(), integer(), integer(), 1, 1)
  }, "data must be finite")
})

test_that("error for Inf data", {
  expect_error({
    LOPART::LOPART_interface(Inf, integer(), integer(), integer(), 1, 1)
  }, "data must be finite")
})

test_that("error for -Inf data", {
  expect_error({
    LOPART::LOPART_interface(-Inf, integer(), integer(), integer(), 1, 1)
  }, "data must be finite")
})
