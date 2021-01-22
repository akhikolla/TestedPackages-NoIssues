context("plots")

test_that("test partial plot", {
   skip("blabla")
   
   Q <- sin(1:10)
   C <- c(rnorm(50), cos(1:10))
   
   tmp <- IncDTW::dtw(Q = rev(Q), C = rev(C), return_QC = TRUE, return_wp = TRUE)
   par <- dtw_partial(tmp, partial_Q = FALSE, partial_C = TRUE, reverse = TRUE)
   tmp$dm <- tmp$dm[nrow(tmp$dm):1, ncol(tmp$dm):1]
   # causes error because NA in dm causes error in BACKTRACK
   plot(tmp, Q = Q, C = C, partial = par)
   
   
   partial <- par
   new_wp <- BACKTRACK_cpp(tmp$dm[  partial$rangeQ[1]:partial$rangeQ[2],
                                  partial$rangeC[1]:partial$rangeC[2]  ])
   rev(new_wp$ii)
   rev(new_wp$jj)
})

