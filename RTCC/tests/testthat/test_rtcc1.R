context("Testing rtcc1")

test_that("rtcc1", {
  data(metadata)
  data(group_information)
  data(group_information)
  RNGversion("3.5")
  set.seed(999999)
  expect_equal(sum(rtcc1(group_information, table_presence_absence, metadata, 2:11, 10)), 306.7)
})
