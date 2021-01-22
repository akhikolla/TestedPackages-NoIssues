context("Testing rtcc2")

test_that("rtcc2", {
  data(metadata)
  data(group_information)
  data(table_presence_absence)
  RNGversion("3.5")
  set.seed(999999)
  expect_equal(round(sum(rtcc2(group_information, table_presence_absence, metadata, group_information$sums, 9, 12, 13, 2, 2, 1, model = 1)), digits = 4), 973.4132)
})
