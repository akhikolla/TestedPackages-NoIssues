context("test running dtw with zscale")

goal <- function(h, x, ws, dm, sp, scale_type){
   nx <- length(x)
   nh <- length(h)
   hscale <- IncDTW::scale(h, type = scale_type)
   sapply(1:(nx-nh+1), function(i){
      y <- IncDTW::scale(x[i:(i+nh-1)], type = scale_type)
      IncDTW::dtw2vec(y, hscale, dist_method = dm,
                      step_pattern = sp, ws = ws)$dist
   })
}


noise <- function(i) cumsum(rnorm(i))


test_that("equal sapply univariate scale1, scale z", {
   
   dm <- "norm1"
   sp <- "symmetric1"
   WS <- 10
   
   h <- noise(10)
   hscale <- IncDTW::scale(h, type="z")
   x <- c(noise(10), h, noise(10))
   
   ist <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "z", threshold = NULL, lower_bound = F)
   soll <- goal(h, x, WS, dm, sp, "z")
   
   
   # sapply(1:(length(x)-length(h)+1), function(i){mean(x[i:(i+length(h)-1)])})
   # sapply(1:(length(x)-length(h)+1), function(i){sd(x[i:(i+length(h)-1)]) })
   # sapply(1:(length(x)-length(h)+1), function(i){var(x[i:(i+length(h)-1)]) })
   # # sapply(1:(length(x)-length(h)+1), function(i){sd(x[i:(i+length(h)-1)])*length(h)/(length(h)-1) })
   # 
   # y <- x[1:(1+length(h)-1)]
   # var(y)
   # myvar0 <- function(x){
   #    mu <- mean(x)
   #    ((1/(length(x)-1)) * sum((x - mu)^2))
   # }
   # myvar2 <- function(x){
   #    mu2 <- mean(x)^2
   #    ss <- sum(x^2 - mu2)
   #    ss/(length(x)-1)
   # }
   # myvar3 <- function(x){
   #    N <- length(x)
   #    mu2 <- mean(x)^2
   #    ss <- sum(x^2)
   #    ss/(N-1) - mu2*N/(N-1) 
   # }
   # c(myvar0(y), myvar2(y), myvar3(y), var(y))
   
   expect_equal(ist$dist, soll)
   
})





test_that("lot-mode", {
   
   dm <- "norm1"
   sp <- "symmetric1"
   WS <- 10
   
   h <- noise(10)
   hscale <- IncDTW::scale(h, type="z")
   x1 <- c(noise(10), h, noise(10))
   x2 <- c(noise(10), h, noise(10))
   x3 <- c(noise(10), h, noise(10))
   x <- c(x1, Inf, x2, Inf, x3)
   
   
   ist <- rundtw(Q = h, C = list(x1, x2, x3), dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "01", threshold = NULL, lower_bound = F)
   ist <- unlist(ist$dist)
   soll <- goal(h, x, WS, dm, sp, "01")
   soll <- soll[!is.na(soll)]
   expect_equal(ist, soll)
   
   
   
   ist <- rundtw(Q = h, C = list(x1, x2, x3), dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "z", threshold = NULL, lower_bound = F)
   ist <- unlist(ist$dist)
   soll <- goal(h, x, WS, dm, sp, "z")
   soll <- soll[!is.na(soll)]
   expect_equal(ist, soll)
   
})





test_that("lot-mode kNN", {
   
   dm <- "norm1"
   sp <- "symmetric1"
   WS <- 10
   
   h <- noise(10)
   hscale <- IncDTW::scale(h, type="z")
   x1 <- c(noise(20))
   x2 <- c(noise(10), h, noise(10))
   x3 <- c(noise(22))
   x4 <- c(noise(10), h, noise(10))
   x <- c(x1, Inf, x2, Inf, x3, Inf, x4)
   
   ret <- rundtw(Q = h, C = list(x1, x2, x3, x4), dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "z", threshold = NULL, lower_bound = F, k = 3)
   
   counter <- 0
   for(i in 1:length(ret$knn_indices)){
      if(ret$dist[[ ret$knn_list_indices[i] ]][ret$knn_indices[i]] == ret$knn_values[i]){
         counter <- counter + 1
      }
   }
   ist <- counter
   soll <- length(ret$knn_indices)
   
   
   expect_equal(ist, soll)
   
})
