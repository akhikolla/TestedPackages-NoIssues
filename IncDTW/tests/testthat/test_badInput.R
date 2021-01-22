context("what happens with bad input")

test_that("test initial_dim_check()", {
   #skip("blabla")
   Q <- cumsum(rnorm(100))
   C <- matrix(rnorm(240), ncol=3, nrow=80)
   expect_error( dtw2vec(Q = Q, C = C))
   
   Q <- matrix(rnorm(240), ncol=2)
   C <- matrix(rnorm(240), ncol=3, nrow=80)
   expect_error( dtw2vec(Q = Q, C = C))
   
   
   Q <- cumsum(rnorm(100))
   C <- "cm"
   expect_error( dtw(Q = Q, C = C))
   
   Q <- matrix(rnorm(100,10,10))
   C <- "cdrm"
   expect_error( dtw(Q = Q, C = C))
   
})

