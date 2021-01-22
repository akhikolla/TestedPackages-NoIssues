context("dba")
# library(dtw)

test_that("correct length of output", {
   WS <- 15
   lot <- lapply(1:10, function(i){
      matrix(rnorm(sample(seq(20, 30, 2), 1), log(i+1)), ncol=2)
   })
   m0 <- lot[[1]]
   ITERMAX <- 10
   tmp <- dba(m0 = m0, lot = lot, iterMax = ITERMAX, eps = NULL, dist_method = "norm2", ws = WS,
                iter_dist_method = "max")
   m1 <- tmp$iterations[[ ITERMAX + 1 ]]
   
   sum_dist_0 <- sum(sapply(lot, function(ll){
      dtw2vec(ll, m0, dist_method = "norm2", ws = WS)$distance
   }))
   
   sum_dist_1 <- sum(sapply(lot, function(ll){
      dtw2vec(ll, m1, dist_method = "norm2", ws = WS)$distance
   }))
   
   expect_equal(length(tmp$iterDist_m2m), ITERMAX )
   expect_equal(length(tmp$iterations), ITERMAX + 1)
   expect_equal(sum_dist_0 >= sum_dist_1, TRUE)
   
})

