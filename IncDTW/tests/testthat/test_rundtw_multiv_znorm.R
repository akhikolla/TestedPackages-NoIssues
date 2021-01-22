context("test running dtw with zscale")

goal <- function(h, x, ws, dm, sp){
   nx <- nrow(x)
   nh <- nrow(h)
   hscale <- IncDTW::scale(h, type = "z")
   sapply(1:(nx-nh+1), function(i){
      y <- IncDTW::scale(x[i:(i+nh-1), , drop = F], type = "z")
      IncDTW::dtw2vec(y, hscale, dist_method = dm,
                      step_pattern = sp, ws = ws)$dist
   })
}


noise <- function(i, nc) matrix(cumsum(rnorm(i*nc)), nrow =i, ncol=nc)

test_that("norm1_sym1", {
   
   dm <- "norm1"
   sp <- "symmetric1"
   nc <- 3
   WS <- 10
   
   h <- noise(10, nc)
   x <- rbind(noise(10, nc), h, noise(10, nc))
   
   ist <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "z", threshold = NULL, lower_bound = F)
   soll <- goal(h, x, WS, dm, sp)
   
   expect_equal(ist$dist, soll)
})


test_that("norm1_sym1 kNN", {
   
   dm <- "norm1"
   sp <- "symmetric1"
   nc <- 3
   WS <- 10
   
   h <- noise(10, nc)
   x <- rbind(noise(10, nc), h, noise(10, nc), h, noise(10, nc), h, noise(10, nc), h, noise(10, nc))
   
   ist <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "z", k = 10)
   soll <- goal(h, x, WS, dm, sp)
   ix_nan <- which(!is.na(ist$dist))
   expect_equal(ist$dist[ix_nan], soll[ix_nan])
})



test_that("lot_mode", {
   
   dm <- "norm1"
   sp <- "symmetric1"
   nc <- 3
   WS <- 10
   
   h <- noise(10, nc)
   x1 <- rbind(noise(10, nc), h, noise(10, nc))
   x2 <- rbind(noise(10, nc), h, noise(10, nc))
   x3 <- rbind(noise(10, nc), h, noise(10, nc))
   x4 <- rbind(noise(10, nc), h, noise(10, nc))
   
   i0 <- matrix(Inf, ncol=3)
   x <- rbind(x1, i0, x2, i0, x3, i0, x4)
   
   ist <- rundtw(Q = h, C = list(x1, x2, x3, x4), dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "z", threshold = NULL, lower_bound = FALSE)
   soll <- goal(h, x, WS, dm, sp)
   soll <- soll[!is.na(soll)]
   ist <- unlist(ist$dist)
   expect_equal(ist, soll)
})



test_that("norm2_sym2", {
   
   dm <- "norm2"
   sp <- "symmetric2"
   nc <- 3
   WS <- 10
   
   h <- noise(10, nc)
   x <- rbind(noise(10, nc), h, noise(10, nc))
   
   ist <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "z", threshold = NULL, lower_bound = F)
   soll <- goal(h, x, WS, dm, sp)
   
   expect_equal(ist$dist, soll)
})




