context("dtw_disvec")
# library(dtw)

test_that("correct values univariate", {
   
   WS <- 11
   eps <- 1/10^7
   lot <- lapply(1:10, function(i){
      rnorm(sample(seq(20, 30, 2), 1), log(i+1))
   })
   Q <- rnorm(20)
   dm1 <- dtw_disvec(Q, lot, ws = WS, normalize = TRUE, dist_method = "norm1", ncores=1)
   dm2 <- dtw_disvec(Q, lot, ws = WS, normalize = TRUE, dist_method = "norm1", ncores=2)
   
   byhand <- sapply(1:10, function(i){
      dtw2vec(Q = Q, lot[[ i ]], dist_method = "norm1", ws = WS)$normalized_distance
   })
   sum(abs(byhand-dm1$disvec))
   sum(abs(byhand-dm2$disvec))
   
   expect_equal(sum(abs(byhand-dm1$disvec)) < eps, TRUE)
   expect_equal(sum(abs(byhand-dm2$disvec)) < eps, TRUE)
})



test_that("correct values multivariate", {
   
   WS <- 6
   eps <- 1/10^7
   lot <- lapply(1:10, function(i){
      matrix(rnorm(sample(seq(20, 30, 2), 1), log(i+1)), ncol = 2)
   })
   Q <- matrix(rnorm(20), ncol = 2)
   dm1 <- dtw_disvec(Q, lot, ws = WS, normalize = TRUE, dist_method = "norm1", ncores=1)
   dm2 <- dtw_disvec(Q, lot, ws = WS, normalize = TRUE, dist_method = "norm1", ncores=2)
   
   
   byhand <- sapply(1:10, function(i){
      dtw2vec(Q = Q, lot[[ i ]], dist_method = "norm1", ws = WS)$normalized_distance
   })
   # sum(abs(byhand-dm1$disvec))
   # sum(abs(byhand-dm2$disvec))
   
   expect_equal(sum(abs(byhand-dm1$disvec)) < eps, TRUE)
   expect_equal(sum(abs(byhand-dm2$disvec))< eps, TRUE)
})


test_that("pass names of input", {
   lot <- lapply(1:4, function(i){
      cumsum(rnorm(10, i))
   })
   names(lot) <- letters[1:4]
   x <- cumsum(rnorm(10))
   
   ret <- dtw_disvec(Q=x, lot = lot,  dist_method="norm2")
   
   expect_equal(is.null(labels(ret$disvec)), FALSE)
   
})