context("Equal dtw package")
library(dtw)


test_that("Bub: dtw normalized mv", {
   Q <- matrix(rnorm(100), ncol=2)
   C <- matrix(rnorm(800), ncol=2)
   
   d_norm1 <- function(x,y){
      sum(abs(x-y))
   }
   cm0 <- IncDTW::cm(Q, C, dist_method = "norm1")
   cm1 <- IncDTW::cm(Q, C, dist_method = d_norm1)
   expect_equal(cm1, cm0)
   
   
   d_norm2 <- function(x,y){
      sqrt(sum((x-y)^2))
   }
   cm0 <- IncDTW::cm(Q, C, dist_method = "norm2")
   cm1 <- IncDTW::cm(Q, C, dist_method = d_norm2)
   expect_equal(cm1, cm0)
   
   
   d_norm22 <- function(x,y){
      sum((x-y)^2)
   }
   cm0 <- IncDTW::cm(Q, C, dist_method = "norm2_square")
   cm1 <- IncDTW::cm(Q, C, dist_method = d_norm22)
   expect_equal(cm1, cm0)
   
})


test_that("custom dist function",{
   rw <- function(nn) sample(1:64, nn, replace=T)
   
   h <- rw(100)
   wc <- log2(3)
   card <- 64
   
   
   x <- rw(100)
   logp <- -c(63:1, 0, 1:63)
   
   
   
   d_mdl <- function(x, y, log_prob){
      as.numeric(-log_prob[x-y + 64])
   }
   
   d_mdl(x[1], h[1], logp)
   
   tmp <- IncDTW::cm(x, h, dist_method = d_mdl, log_prob = logp)
   counter <- 0
   for(i in 1:nrow(tmp)){
      for(j in 1:ncol(tmp)){
         if(tmp[i,j] == d_mdl(x[i], h[j], logp)) counter <- counter + 1
      }
   }
   
   expect_equal(counter, prod(dim(tmp)))
})




