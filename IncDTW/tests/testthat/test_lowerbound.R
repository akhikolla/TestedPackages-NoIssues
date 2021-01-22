context("lowerbound")


test_that("lowerbound smaller than dtw distance, univariate", {
   
   
   counter <- 0
   for(i in 1:30){
      
      ws <- sample(2:40, size = 1)
      # dist_method <- sample(c("norm1", "norm2", "norm2_square"), size = 1)
      dist_method <- "norm1"
      nh <- 50
      nx <- 50
      h <- cumsum(rnorm(nh))
      x <- cumsum(rnorm(nx))
      h.z <- IncDTW::scale(h, "z")
      x.z <- IncDTW::scale(x, "z")
      
      lb.z <- lowerbound(C = x.z, ws = ws, scale ="none", dist_method = dist_method, Q = h.z)
      lb <- lowerbound(C = x, ws = ws, scale ="z", dist_method = dist_method, Q = h)
      d1 <- dtw2vec(Q = h.z, C = x.z, step_pattern = "symmetric1", dist_method = dist_method, ws = ws)$distance
      d2 <- dtw2vec(Q = h.z, C = x.z, step_pattern = "symmetric2", dist_method = dist_method, ws = ws)$distance
      
      if(lb.z > d1) counter <- counter + 1
      if(lb.z > d2) counter <- counter + 1
      if(lb != lb.z) counter <- counter + 1
      # expect_gt(d, lb.z)
      # expect_equal(lb, lb.z)
   }
   expect_equal(counter, 0)
   
})




test_that("lowerbound smaller than dtw distance, multivariate", {
   
   
   counter <- 0
   for(i in 1:30){
      
      ws <- sample(2:40, size = 1)
      dist_method <- sample(c("norm1", "norm2", "norm2_square"), size = 1)
      # dist_method <- "norm1"
      nh <- 50
      nx <- 50
      nc <- 3
      h <- matrix(cumsum(rnorm(nh * nc)), ncol = nc)
      x <- matrix(cumsum(rnorm(nx * nc)), ncol = nc)
      h.z <- IncDTW::scale(h, "z")
      x.z <- IncDTW::scale(x, "z")
      
      lb.z <- lowerbound(C = x.z, ws = ws, scale ="none", dist_method = dist_method, Q = h.z)
      lb <- lowerbound(C = x, ws = ws, scale ="z", dist_method = dist_method, Q = h)
      d1 <- dtw2vec(Q = h.z, C = x.z, step_pattern = "symmetric1", dist_method = dist_method, ws = ws)$distance
      d2 <- dtw2vec(Q = h.z, C = x.z, step_pattern = "symmetric2", dist_method = dist_method, ws = ws)$distance
      
      if(lb.z > d1) counter <- counter + 1
      if(lb.z > d2) counter <- counter + 1
      if(lb != lb.z) counter <- counter + 1
      # expect_gt(d, lb.z)
      # expect_equal(lb, lb.z)
   }
   expect_equal(counter, 0)
   
})


test_that("lowerbound smaller than dtw distance - test", {
   # test it with the function cpp_get_tube
   # which is not exported to an R fucntion but also 
   # available after devtools::load_all(".")
   # => the final counter should be 0!
  skip("")
   
   
   
   ws <- 10
   nh <- 100
   nx <- 100
   h <- cumsum(rnorm(nh))
   x <- cumsum(rnorm(nx))
   tube <- cpp_get_tube(h, ws);
   j0 <- 0
   jsup <- 0+nh
   get_lb(tube, x, 0, jsup)
   dtw2vec(h, x[(j0+1):(jsup)], ws=ws, dist_method = "norm1", step_pattern = "symmetric1")$distance
   
   
   
   plot(h, type="l")
   lines(tube[,1], col="blue")
   lines(tube[,2], col="red")
   
   
   
   
   counter <- 0
   for(i in 1:1000){
      myseed <- sample(10^4, size = 1)
      set.seed(myseed)
      
      nh <- 100
      nx <- 100
      ws <- sample(nh, size = 1)
      h <- cumsum(rnorm(nh))
      x <- cumsum(rnorm(nx))
      tube <- cpp_get_tube(h, ws);
      j0 <- 0
      jsup <- 0+nh
      lb <- get_lb(tube, x, 0, jsup)
      d <- dtw2vec(h, x[(j0+1):(jsup)], ws=ws, dist_method = "norm1", step_pattern = "symmetric1")$distance
      
      if(lb > d){
         print(myseed)
         counter <- counter +1
      }
   }
   print(counter)
   
   
})

