context("test running dtw")

goal <- function(h, x, ws, dm, sp){
   nx <- nrow(x)
   nh <- nrow(h)
   hscale <- IncDTW::scale(h, type = "01")
   sapply(1:(nx-nh+1), function(i){
      y <- IncDTW::scale(x[i:(i+nh-1), , drop = F], type = "01")
      IncDTW::dtw2vec(y, hscale, dist_method = dm,
                      step_pattern = sp, ws = ws)$dist
   })
}

goal_noscale <- function(h, x, ws, dm, sp){
   nx <- nrow(x)
   nh <- nrow(h)
   sapply(1:(nx-nh+1), function(i){
      y <- x[i:(i+nh-1), , drop = F]
      IncDTW::dtw2vec(y, h, dist_method = dm,
                      step_pattern = sp, ws = ws)$dist
   })
}

noise <- function(i, nc) matrix(cumsum(rnorm(i*nc)), nrow =i, ncol=nc)

test_that("correct multivariate norm1", {
  
   dm <- "norm1"
   sp <- "symmetric1"
   nc <- 3
   WS <- 10
   
   h <- noise(10, nc)
   x <- rbind(noise(10, nc), h, noise(10, nc))
   
   # DO scale
   ist <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "01", threshold = NULL, lower_bound = F)
   soll <- goal(h, x, WS, dm, sp)
   
   expect_equal(ist$dist, soll)
   
   # DONT scale
   ist <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "none", threshold = NULL, lower_bound = F)
   soll <- goal_noscale(h, x, WS, dm, sp)
   
   expect_equal(ist$dist, soll)
   
})


test_that("correct multivariate norm2", {
   
   dm <- "norm2"
   sp <- "symmetric1"
   nc <- 3
   WS <- NULL

   h <- noise(10, nc)
   hscale <- IncDTW::scale(h, type="01")
   x <- rbind(noise(10, nc), h, noise(10, nc))
   
   # DO scale
   ist <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "01", threshold = NULL, lower_bound = F)
   soll <- goal(h, x, WS, dm, sp)
   
   expect_equal(ist$dist, soll)
   
   # DONT scale
   ist <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "none", threshold = NULL, lower_bound = F)
   soll <- goal_noscale(h, x, WS, dm, sp)
   
   expect_equal(ist$dist, soll)
   
   
})


test_that("correct multivariate sym2", {
   
   dm <- "norm2"
   sp <- "symmetric2"
   nc <- 3
   WS <- 100
   
   h <- noise(10, nc)
   hscale <- IncDTW::scale(h, type="01")
   x <- rbind(noise(10, nc), h, noise(10, nc))
   
   # DO scale
   ist <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "01", threshold = NULL, lower_bound = F)
   soll <- goal(h, x, WS, dm, sp)
   
   expect_equal(ist$dist, soll)
   
   # DONT scale
   ist <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "none", threshold = NULL, lower_bound = F)
   soll <- goal_noscale(h, x, WS, dm, sp)
   
   expect_equal(ist$dist, soll)
   
})


test_that("expected result, multivariate", {
   
   dm <- "norm1"
   sp <- "symmetric1"
   nc <- 3
   WS <- NULL
   deform <- function(x) (x + rnorm(1, 0, 100)) * abs(rnorm(1, 0, 100))
   
   h <- noise(10, nc)
   x <- rbind(noise(10, nc), h, noise(10, nc), deform(h), noise(10, nc), deform(h))
   soll <- c(11, 31, 51)
   
   # DONT lowerbound
   ist <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, scale = "01",
                 threshold = NULL, lower_bound = FALSE)
   ist <- which(ist$dist < 10^(-10))
   
   expect_equal(ist, soll)
   
   
   # DO lowerbound
   ist <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, scale = "01",
                 threshold = NULL, lower_bound = TRUE)
   ist <- which(ist$dist < 10^(-10))
   
   expect_equal(ist, soll)
   
})



test_that("kNN rundtw", {
   dm <- "norm1"
   sp <- "symmetric1"
   WS <- NULL
   nc <- 3
   goal_knn <- function(hscale, x, sp, ws = 10, kNNk){
      allv <- rundtw(h, x, dist_method = dm, step_pattern = sp, ws = 10, 
                     threshold = NULL, k = 0, lower_bound = FALSE, scale = "01")
      dd <- allv$dist
      best_i <- integer()
      for(i in 1:kNNk){
         ix <- which.min(dd)
         if(dd[ix] < Inf){
            best_i <- c(best_i, ix)
            dd[max(1, (1 + ix - nh)):min((ix + nh - 1), length(dd))] <- Inf   
         }
      }
      best_i
   }
   
   
   nx <- 500
   nh <- 20
   nfits <- 10
   
   nn <- nx - nfits * nh# nnoise
   nn <- nn/nfits
   
   h <- noise(nh, nc)
   hscale <- IncDTW::scale( h , type="01")
   x <- noise(0,nc)
   for(i in 1:nfits){
      x <- rbind(x, noise(nn, nc), h)
   }
   
   # rundtw(hscale, x,dm,  sp, k = 5, ws = 10, threshold = Inf,1,1,1)# DO scale
   # cpp_rundtw(hscale, x, sp, ws = 10, threshold = Inf, kNNk = 10, do_norm = 1, use_ea = 1, use_lb = 1)
   
   k <- nfits + 2
   
   
   tmp <- rundtw(h, x, dist_method = dm, step_pattern = sp, ws = 10, threshold = Inf, k = k,
                 lower_bound = TRUE, scale = "01")
   ist <- sort(tmp$knn_indices)
   
   soll <- sort( goal_knn(hscale, x, sp, ws=10, kNNk = k) )
   
   
   expect_equal(ist, soll)
})





test_that("run time multiv", {
   skip("runtime comparison")
   
   dm <- "norm1"
   sp <- "symmetric1"
   WS <- NULL
   noise <- function(nr, nc) matrix(cumsum(rnorm(nr*nc)), nrow =nr, ncol=nc)
   # deform <- function(x) (x + rnorm(1, 0, 100)) * abs(rnorm(1, 0, 100))
   
   # speed comparison
   foo_2vec <- function(hscale, x, ws){
      nh <- nrow(hscale)
      nx <- nrow(x)
      sapply(1:(nx-nh+1), function(i){
         y <- IncDTW::scale(x[i:(i+nh-1), , drop=F], type = "01")
         IncDTW::dtw2vec(y, hscale, dist_method = dm,
                         step_pattern = sp, ws = ws)$dist
      })
   }
   
   foo_cm <- function(hscale, x, ws){
      nh <- nrow(hscale)
      nx <- nrow(x)
      sapply(1:(nx-nh+1), function(i){
         y <- IncDTW::scale(x[i:(i+nh-1), , drop=F], type = "01")
         cm_tmp <- IncDTW::cm(y, hscale, dist_method = dm, ws = ws)
         IncDTW::dtw2vec(cm_tmp, C = "cm", step_pattern = sp, ws = ws)$dist
      })
   }
   
   
   NH <- c(0.2, 0.3)
   NX <- c(100, 500, 750)
   
   # NH <- seq(0.1, 0.5, 0.05)
   # NX <- c(100, 200, 300, 500, 1000)
   NC <- c(2)
   
   df <- data.frame(expr = character(),
                    min = numeric(), 
                    mean = numeric(), 
                    median = numeric(), 
                    max = numeric(), 
                    neval = integer(), 
                    rel = numeric(), 
                    nh = integer(),
                    nx = integer(),
                    nfits = integer(),
                    nc = integer()
   )
   
   data.table::fwrite(df, file = "comp_time_test_multiv.txt")
   tic <- Sys.time()
   for(nhh in NH){
      print(nhh); print(Sys.time())
      
      for(nx in NX){
         print(c(nhh, nx)); print(Sys.time())
         
         for(nc in NC){
            print(c(nhh, nx, nc)); print(Sys.time())
            
            nh <- nhh*nx
            nfits0 <- 10
            nfits1 <-  min(nfits0, floor(nx/nh))
            
            nfits <- sample(nfits1, size = 1)
            
            nn <- nx - nfits * nh# nnoise
            nn <- nn/nfits
            
            h <- noise(nh, nc)
            hscale <- IncDTW::scale( h , type="01")
            x <- matrix(numeric(), ncol=nc)
            for(i in 1:nfits){
               x <- rbind(x, noise(nn, nc), h)
            }
            
            
            
            mic <- microbenchmark::microbenchmark(foo_cm(hscale, x, WS),
                                                  foo_2vec(hscale, x, WS),
                                                  rundtw(hscale, x, dm, sp, "01", WS, NULL),
                                                  rundtw(hscale, x, dm, sp, "01", WS, Inf),
                                                  rundtw(hscale, x, dm, sp, "none", WS, NULL),
                                                  rundtw(hscale, x, dm, sp, "none", WS, Inf),
                                                  times = 10)
            df0 <- maxaR::rel_microbenchmark(mic, 
                                             cols = c("expr", "min", "mean", "median", "max", "neval"), 
                                             rel_expr = "rundtw(hscale, x, dm, sp, T, WS, Inf)")$df
            df0$nh <- nh
            df0$nx <- nx
            df0$nfits <- nfits
            df0$nc <- nc
            
            data.table::fwrite(df0, file = "comp_time_test_multiv.txt", append = T, row.names = F, col.names = F)
            df <- rbind(df, df0)
         }
      }
   }
   print(Sys.time())
   
   require(ggplot2)
   
   df <- data.table::fread(file = "comp_time_test_multiv.txt")
   ggplot(df) + geom_line(aes(x = nh/nx, y = rel, group = expr, col = expr), size = 2)+
      facet_grid(~nx) + scale_y_log10()
   
   ggplot(df) + geom_line(aes(x = nh/nx, y = median, group = expr, col = expr), size = 2)+
      facet_grid(~nx) + scale_y_log10()
   
   df[, list(mean(median)), by = c( "nfits")]
   
})