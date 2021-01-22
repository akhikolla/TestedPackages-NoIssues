context("planedtw and friends")


test_that("in-/decrement", {
   
   rw <- function(nn) cumsum(rnorm(nn))
   Q <- sin(1:100)
   C <- Q[1:90] + rnorm(90, 0, 0.1)
   WS <- 40
   # initial calculation
   x <- initialize_plane(Q, C, ws = WS)
   
   # incremental calculation for new observations
   y1 <- Q[91:92] + rnorm(2, 0, 0.1)# new observations
   x <- increment(x, newObs = y1)
   y2 <- Q[93:95] + rnorm(3, 0, 0.1)# new observations
   x <- increment(x, newObs = y2)
   y3 <- c(Q[96:100] + rnorm(5, 0, 0.1), rw(10))# new observations
   x <- increment(x, newObs = y3)
   
   ist <- x$control$nC
   soll <- 110
   expect_equal(ist, soll)
   
   ist <- x$normalized_distance 
   soll <- dtw2vec(Q, c(C, y1, y2, y3) , ws = WS)$normalized_distance
   expect_equal(ist, soll)
   
   x <- decrement(x, direction = "C", refresh_dtw = TRUE)
   ist <- x$normalized_distance 
   soll <- dtw2vec(Q, c(C, y1, y2, y3[1:5]) , ws = WS)$normalized_distance
   expect_lte(ist, soll) # 'ist' can be smaller, if the partial found fit is of diff length
   
})


test_that("reverse increment", {
   
   
   #--- 2. example: first increment, then reverse increment
   rw <- function(nn) cumsum(rnorm(nn))
   Q <- rw(100)
   C <- Q[11:90] + rnorm(80, 0, 0.1)
   WS <- 40
   # initial calculation
   x <- initialize_plane(Q, C, ws = WS)
   
   # incremental calculation for new observations
   y1 <- Q[91:100] + rnorm(10, 0, 0.1)# new observations
   x <- increment(x, newObs = y1)
   x <- reverse(x)
   y2 <- Q[10:6] + rnorm(5, 0, 0.1)# new observations in reverse order
   x <- increment(x, newObs = y2)
   y3 <- Q[5:1] + rnorm(5, 0, 0.1)# again new obs in reverse order
   x <- increment(x, newObs = y3)
   #x <- reverse(x)
   ist <- x$distance 
   soll <- dtw2vec(rev(Q), rev(c(rev(y3), rev(y2), C, y1)),ws = WS)$distance
   
   expect_equal(ist, soll)
   
   
})


test_that("refresh", {
   
   
   rw <- function(nn) cumsum(rnorm(nn))
   Q <- rw(100)
   C <- Q[11:90] + rnorm(80, 0, 0.1)
   WS <- 40
   # initial calculation
   x <- initialize_plane(Q, C, ws = WS)
   x <- refresh(x)
   
   ist <- x$normalized_distance
   soll <- dtw2vec(Q, C, ws = WS)$normalized_distance
   expect_equal(ist, soll)
   
   
   x <- decrement(x, refresh_dtw = FALSE)
   x <- refresh(x)
   ist <- x$normalized_distance
   soll <- dtw2vec(Q, C[1:x$control$nC], ws = WS)$normalized_distance
   expect_equal(ist, soll)
   
})



test_that("primitive decrement", {
   
   
   rw <- function(nn) cumsum(rnorm(nn))
   Q <- rw(100)
   C <- Q[11:90] + rnorm(80, 0, 0.1)
   WS <- 100
   # initial calculation
   x <- initialize_plane(Q, C, ws = WS)
   x <- decrement(x, nC = 20)
   expect_equal(x$gcm_lc_new, NULL)
   expect_equal(x$gcm_lr_new, NULL)
   expect_equal(x$control$nC, 20)
   
   x <- initialize_plane(Q, C, ws = WS)
   x <- decrement(x, nC = 20, refresh_dtw = TRUE)
   
   expect_equal(x$normalized_distance, dtw2vec(Q, C[1:20], ws = WS)$normalized_distance)
   expect_equal(x$control$nC, 20)
   
})



test_that("partial decrement", {
   
   rw <- function(nn) cumsum(rnorm(nn))
   Q <- rw(100)
   C <- c(Q[1:90] + rnorm(90, 0, 0.1), rnorm(20))
   WS <- 30
   
   # initial calculation
   x <- initialize_plane(Q, C, ws = WS)
   par <- dtw_partial(x, partial_Q = FALSE, partial_C = TRUE)
   x1 <- decrement(x, refresh_dtw = TRUE)
   x2 <- decrement(x, nC = par$rangeC[2], refresh_dtw = TRUE)
   expect_equal(x1, x2)
   
   x1 <- decrement(x, refresh_dtw = FALSE)
   x2 <- decrement(x, nC = par$rangeC[2], refresh_dtw = FALSE)
   expect_equal(x1, x2)
   
   x <- decrement(x, direction = "C")
   expect_equal(x$gcm_lc_new, NULL)
   expect_equal(x$gcm_lr_new, NULL)
   
   x <- initialize_plane(Q, C, ws = WS)
   x <- decrement(x, nC = 75, refresh_dtw = TRUE)
   
   tmp <- dtw2vec(Q, C[1:75], ws = WS)
   expect_equal(x$normalized_distance, tmp$normalized_distance)
   expect_equal(x$distance, tmp$distance)
   expect_equal(x$control$nC, 75)
   

})
