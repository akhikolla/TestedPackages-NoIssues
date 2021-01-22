context("Testing dist1.cpp")

test_that("dist1", {
  expect_equal(dist1(1:10), sum(stats::dist(1:10)))
})


