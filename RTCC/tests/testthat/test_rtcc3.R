context("Testing rtcc3")

test_that("rtcc3", {
  data(metadata)
  data(group_information)
  data(group_information)
  RNGversion("3.5")
  set.seed(999999)
  expect_equal(sum(rtcc3(group_information, table_presence_absence, metadata, group_information$sums, 9, 12, 13, 2, 1, 1, model = 1)), 879.784441)
})
