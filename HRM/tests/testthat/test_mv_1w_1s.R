context("1 Wholeplot, 1 Subplot")
true_fH <- c(2.81416895836941 ,  2.81416895836941 ,  2.81416895836941 , 11.2566758334776, 1 ,  1 ,  1 , 36,  2.90156459950675 ,  2.90156459950675 ,  2.90156459950675, 104.456325582243)
true_fG <- c(97.5996024122159 ,  97.5996024122159 ,  97.5996024122159 , 241.15186,  110.254489338651 ,  110.254489338651 ,  110.254489338651 , 75.25449, 110.254489338651 ,  110.254489338651 ,  110.254489338651, 218.80741 )
true_test <- c(0.636499758309033 ,  0.645610906141085 ,  0.6349818739943 , 0.5212962, 3353.90095426749 ,  3369.07643764296 ,  3357.32491154643 , 13984.51, 3.63813296984172 ,  3.68569441019603 ,  3.65440661642982, 1.322141)
result <- HRM::hrm_test(value ~ group*region, subject = subject, variable = variable, data = EEG)$result

test_that("function hrm_test MVHRM", {
  # expect_equal(result$fH, true_fH, tol = 1e-4)
  # expect_equal(result$fG, true_fG, tol = 1e-4)
  # expect_equal(result$test, true_test, tol = 1e-4)
  # expect_output(print(HRM::hrm_test(value ~ group*region, subject = subject, variable = variable, data = EEG)))
  # expect_output(summary(HRM::hrm_test(value ~ group*region, subject = subject, variable = variable, data = EEG)))
})
