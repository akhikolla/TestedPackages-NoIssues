context("Multivariate time series")
# library(dtw)

test_that("IncDTW::dtw2vec EQUAL IncDTW::dtw with former cm, norm2", {
   Q <- matrix(rnorm(100),ncol=2)
   C <- matrix(rnorm(100),ncol=2)
   WS <- 5
   WS2 <- 5000
   fooVec <- function(Q,C,ws){ dtw2vec(Q = Q, C = C, dist_method = "norm2", ws = ws)$distance }
   fooMat <- function(Q,C,ws){ 
      cm0 <- IncDTW::cm(Q, C, dist_method = "norm2", ws = ws)
      IncDTW::dtw(Q = cm0, C = "cm", dist_method = "norm2", ws = ws)$distance}
   
   expect_equal(fooVec(Q, C, ws = NULL), fooMat(Q, C, ws = NULL))
   expect_equal(fooVec(Q, C, ws = 0),    fooMat(Q, C, ws = 0))
   expect_equal(fooVec(Q, C, ws = WS),   fooMat(Q, C, ws = WS))
   expect_equal(fooVec(Q, C, ws = WS2),  fooMat(Q, C, ws = WS2))
})



test_that("IncDTW::dtw2vec EQUAL IncDTW::dtw with former cm, norm1", {
   Q <- matrix(rnorm(100),ncol=2)
   C <- matrix(rnorm(100),ncol=2)
   WS <- 5
   WS2 <- 5000
   fooVec <- function(Q,C,ws){ dtw2vec(Q = Q, C = C, dist_method = "norm1", ws = ws)$distance }
   fooMat <- function(Q,C,ws){ 
      cm0 <- IncDTW::cm(Q, C, dist_method = "norm1", ws = ws)
      IncDTW::dtw(Q = cm0, C = "cm", dist_method = "norm1", ws = ws)$distance}
   
   expect_equal(fooVec(Q, C, ws = NULL), fooMat(Q, C, ws = NULL))
   expect_equal(fooVec(Q, C, ws = 0),    fooMat(Q, C, ws = 0))
   expect_equal(fooVec(Q, C, ws = WS),   fooMat(Q, C, ws = WS))
   expect_equal(fooVec(Q, C, ws = WS2),  fooMat(Q, C, ws = WS2))
})



test_that("IncDTW::dtw2vec EQUAL dtw::dtw __ multivariate, norm2", {
   dist_method <- "norm2"
   dist.method <- "Euclidean"
   
   #--- univariate
   
   Q <- rnorm(100)
   C <- rnorm(100)
   ##--- no ws
   d0 <- dtw::dtw(Q,C,step.pattern = dtw::symmetric1, 
                  window.type = "sakoe", window.size = Inf, keep.internals = T)
   d1 <- dtw2vec(Q = Q, C = C, ws = NULL, step_pattern = "symmetric1")$distance 
   identical(d0$distance, d1)
   expect_equal(d0$distance, d1)
   
   ##--- with ws
   WS <- 10
   d0 <- dtw::dtw(Q,C,step.pattern = dtw::symmetric2,
                  window.type = "sakoe", window.size = WS, keep.internals = T)
   d1 <- dtw2vec(Q = Q, C = C, step_pattern = "symmetric2", ws = WS)$distance 
   identical(d0$distance, d1)
   expect_equal(d0$distance, d1)
   
   
   
   #--- multivariate
   Qm <- matrix(rnorm(200),ncol=2)
   Cm <- matrix(rnorm(200),ncol=2)
   
   ##--- no ws
   d0 <- dtw::dtw(Qm, Cm, step.pattern = dtw::symmetric1, dist.method = dist.method, 
                  window.type = "sakoe", window.size = Inf, keep.internals = T)
   d1 <- dtw2vec(Q = Qm, C = Cm, step_pattern = "symmetric1", dist_method = dist_method, ws = NULL)$distance 
   identical(d0$distance, d1)
   expect_equal(d0$distance, d1)
   
   ##--- with ws
   WS <- 10
   d0 <- dtw::dtw(Qm, Cm, step.pattern = dtw::symmetric2, dist.method = dist.method, 
                  window.type = "sakoe", window.size = WS, keep.internals = T)
   d1 <- dtw2vec(Q = Qm, C = Cm, step_pattern = "symmetric2", dist_method = dist_method, ws = WS)$distance 
   identical(d0$distance, d1)
   expect_equal(d0$distance, d1)
   
   
})



test_that("IncDTW::dtw2vec EQUAL dtw::dtw __ multivariate, norm1", {
   dist_method <- "norm1"
   dist.method <- "manhattan"
   
   #--- univariate
   
   Q <- rnorm(100)
   C <- rnorm(100)
   ##--- no ws
   d0 <- dtw::dtw(Q,C,step.pattern = dtw::symmetric2, 
                  window.type = "sakoe", window.size = Inf, keep.internals = T)
   d1 <- dtw2vec(Q = Q, C = C, step_pattern = "symmetric2", ws = NULL)$distance 
   identical(d0$distance, d1)
   expect_equal(d0$distance, d1)
   
   ##--- with ws
   WS <- 10
   d0 <- dtw::dtw(Q,C,step.pattern = dtw::symmetric1,
                  window.type = "sakoe", window.size = WS, keep.internals = T)
   d1 <- dtw2vec(Q = Q, C = C, step_pattern = "symmetric1", ws = WS)$distance 
   identical(d0$distance, d1)
   expect_equal(d0$distance, d1)
   
   
   
   #--- multivariate
   Qm <- matrix(rnorm(10),ncol=2)
   Cm <- matrix(rnorm(10),ncol=2)
   
   ##--- no ws
   d0 <- dtw::dtw(Qm, Cm, step.pattern = dtw::symmetric2, dist.method = dist.method, 
                  window.type = "sakoe", window.size = Inf, keep.internals = T)
   d1 <- dtw2vec(Q = Qm, C = Cm, step_pattern = "symmetric2", dist_method = dist_method, ws = NULL)$distance 
   identical(d0$distance, d1)
   expect_equal(d0$distance, d1)
   
   ##--- with ws
   WS <- 100
   d0 <- dtw::dtw(Qm, Cm, step.pattern = dtw::symmetric1, dist.method = dist.method, 
                  window.type = "sakoe", window.size = WS, keep.internals = T)
   d1 <- dtw2vec(Q = Qm, C = Cm, dist_method = dist_method, step_pattern = "symmetric1",
                 ws = WS, threshold = 10^4)$distance 
   identical(d0$distance, d1)
   expect_equal(d0$distance, d1)
   
   
})
