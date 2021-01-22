## Test evaluation functions on real and toy tasks
context("Evaluation functions")
library(wordspace)

qw <- function (x) unlist(strsplit(x, "\\s+", perl=TRUE)) # Perl's qw()


## test small invented multiple-choice task on minimalistic DSM example
Mochi <- data.frame(
  target      = qw("cat   dog    time   time   cause"),
  correct     = qw("dog   cause  animal reason effect"),
  distractor1 = qw("time  reason cat    animal reason"),
  distractor2 = qw("cause cat    dog    cause  time"),
  stringsAsFactors = TRUE # make sure that things don't break with factors
)
Mochi.ok <- c(TRUE, FALSE, TRUE, FALSE, TRUE)
Mochi.best <- qw("dog cat animal animal effect")

test_that("multiple choice task can be evaluated correctly", {
  M <- dsm.score(DSM_TermTerm, score="frequency", transform="log", matrix.only=TRUE)

  Mochi.eval <- eval.multiple.choice(Mochi, M)
  expect_equal(Mochi.eval$TP, sum(Mochi.ok))
  expect_equal(Mochi.eval$FP, sum(!Mochi.ok))
  
  Mochi.details <- eval.multiple.choice(Mochi, M, details=TRUE)
  expect_equal(Mochi.details$correct, Mochi.ok)
  expect_equal(Mochi.details$best.choice, Mochi.best)
})
