context("test running dtw")

goal <- function(h, x, ws, dm, sp){
   nx <- length(x)
   nh <- length(h)
   hscale <- IncDTW::scale(h, type = "01")
   sapply(1:(nx-nh+1), function(i){
      y <- IncDTW::scale(x[i:(i+nh-1)], type = "01")
      IncDTW::dtw2vec(y, hscale, dist_method = dm,
                      step_pattern = sp, ws = ws)$dist
   })
}


noise <- function(i) cumsum(rnorm(i))


test_that("no variance for a sequence", {
   
   dm <- "norm1"
   sp <- "symmetric1"
   WS <- 10
   
   #my_seed <- 53#sample(100, 1); print(my_seed)
   # set.seed(my_seed)
   h <- rep(rnorm(1), 20)
   hscale <- IncDTW::scale(h, type="01")
   x <- c(noise(10), rep(rnorm(1), 20), noise(10))
   
   ret <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "01", threshold = NULL, lower_bound = F)
   ist <- ret$dist[11]
   soll <- 0
   
   expect_equal(ist, soll)
   
   
   
   ret <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "z", threshold = NULL, lower_bound = F)
   ist <- ret$dist[11]
   soll <- 0
   
   expect_equal(ist, soll)
})




test_that("equal sapply univariate norm1", {
  
   dm <- "norm1"
   sp <- "symmetric1"
   WS <- 10
   
   h <- noise(10)
   hscale <- IncDTW::scale(h, type="01")
   x <- c(noise(10), h, noise(10))
   
   ist <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "01", threshold = NULL, lower_bound = F)
   soll <- goal(h, x, WS, dm, sp)
   
   expect_equal(ist$dist, soll)
   
})


test_that("equal sapply univariate sym2", {
   
   dm <- "norm1"
   sp <- "symmetric2"
   WS <- 10
   
   h <- noise(10)
   hscale <- IncDTW::scale(h, type="01")
   x <- c(noise(10), h, noise(10))
   
   ist <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "01", threshold = NULL, lower_bound = F)
   soll <- goal(h, x, WS, dm, sp)
   
   expect_equal(ist$dist, soll)
   
})


test_that("expected result, univariate", {
   
   dm <- "norm1"
   sp <- "symmetric1"
   WS <- NULL
   deform <- function(x) (x + rnorm(1, 0, 100)) * abs(rnorm(1, 0, 100))
   
   h <- noise(10)
   hscale <- IncDTW::scale(h, type="01")
   x <- c(noise(10), h, noise(10), deform(h), noise(10), deform(h))
   soll <- c(11, 31, 51)
   
   # without lower bound
   ist <- rundtw(Q = hscale, C = x, dist_method = dm, step_pattern = sp, 
                 scale = "01", threshold = NULL, lower_bound = FALSE)
   ist <- which(ist$dist < 10^(-10))
   expect_equal(ist, soll)
   
   # with lower bound
   ist <- rundtw(Q = hscale, C = x, dist_method = dm, step_pattern = sp, 
                 scale = "01", threshold = Inf, lower_bound = TRUE)
   ist <- which(ist$dist < 10^(-10))
   
   expect_equal(ist, soll)
   
})


test_that("rundtw and find_peaks", {
   
   noise <- function(nr) cumsum(rnorm(nr))
   
   nx <- 500
   nh <- 30
   nfits <- 5
   
   nn <- nx - nfits * nh# nnoise
   nn <- nn/nfits
   
   h <- noise(nh)
   hscale <- IncDTW::scale( h , type="01")
   x <- noise(0)
   for(i in 1:nfits){
      x <- c(x, noise(nn), h)
   }
   
   # DO scale and DO lower bound
   result <- rundtw(h, x, scale = "01", ws = 10, threshold = Inf, lower_bound = TRUE)
   ## have a look into the result and get best indices with lowest distance
   
   soll_indices <- c(71, 171, 271, 371, 471)
   minima <- find_peaks(result$dist, w = nh)
   ist <- sum(soll_indices %in% minima)
   soll = length(soll_indices)
   
   # plot(result$dist)
   # points(x = minima, y = result$dist[minima], col="red")
   
   
   expect_equal(ist, soll)
})


test_that("kNN rundtw", {
   noise <- function(nr) cumsum(rnorm(nr))
   
   nx <- 500
   nh <- 30
   nfits <- 5
   
   nn <- nx - nfits * nh# nnoise
   nn <- nn/nfits
   
   h <- noise(nh)
   hscale <- IncDTW::scale( h , type="01")
   x <- noise(0)
   for(i in 1:nfits){
      x <- c(x, noise(nn), h)
   }
   
   k <- 10
   # DO scale and DO lower bound
   result <- rundtw(h, x, scale = "01", ws = 10, threshold = Inf, k = k, lower_bound = TRUE)
   ## have a look into the result and get best indices with lowest distance
   
   
   soll <- c(71, 171, 271, 371, 471)
   ord <- order(result$knn_values, decreasing = F)
   ist <- sort(result$knn_indices[ord[1:5]])
   
   # plot(result$dist)
   # points(x = minima, y = result$dist[minima], col="red")
   
   
   expect_equal(ist, soll)
})


test_that("kNN rundtw", {
   
   for(i in 1:1){
   
      myseed <- sample(i + as.numeric(Sys.time())/10^5, size=1)
       # myseed <- 11378
      print(myseed)
      set.seed(myseed)
      nx <- 500
      nh <- sample(20:50, size=1)
      WS <- sample(nh, size=1)
      dm <- "norm1"
      sp <- "symmetric1"
      noise <- function(nr) cumsum(rnorm(nr))
      goal_knn <- function(hscale, x, sp, ws, kNNk){
         allv <- rundtw(h, x, dist_method = dm, step_pattern = sp, ws = ws, threshold = NULL, k = 0,
                        lower_bound = FALSE, scale = "01")
         dd <- allv$dist
         best_i <- integer()
         for(i in 1:kNNk){
            ix <- which.min(dd)
            if(dd[ix] < Inf){
               best_i <- c(best_i, ix)
               dd[max(1, (1 + ix - nh)):min((ix + nh - 1), length(dd))] <- Inf   
            }
         }
         # print(cbind(1:length(dd), allv$dist, dd))
         return(list(indices=best_i, all_values = allv$dist))
      }
      
      
     
      nfits <- sample(10, size=1)
      
      nn <- nx - nfits * nh# nnoise
      nn <- nn/nfits
      
      h <- noise(nh)
      hscale <- IncDTW::scale( h , type="01")
      x <- noise(0)
      for(i in 1:nfits){
         x <- c(x, noise(nn), h)
      }
      
      # rundtw(hscale, x,dm,  sp, k = 5, ws = 10, threshold = Inf,1,1,1)# DO scale
      # cpp_rundtw(hscale, x, sp, ws = 10, threshold = Inf, kNNk = 10, do_scale = 1, use_ea = 1, use_lb = 1)
      
      k <- nfits + 2
      
      
      a <- rundtw(h, x, dist_method = dm, step_pattern = sp, ws = WS, threshold = Inf, k = k,
                    lower_bound = TRUE, scale = "01")
      ist <- sort(a$knn_indices)
      
      b <- goal_knn(hscale, x, sp, ws=WS, kNNk = k)
      soll <- sort( b$indices )
      
      
      # cbind(soll, b$all_values[soll])
      # cbind(ist, b$all_values[ist])
      # 
      # zz <- file("build_ignore/mySink.Rout", open = "wt"); sink(zz); sink(zz, type = "message")
      # tmp <- rundtw(h, x, dist_method = dm, step_pattern = sp, ws = WS, threshold = Inf, k = k,
      #               lower_bound = TRUE, scale = "01", debug = 1)
      # 
      # cpp_kNN_rev(tmp$dist, nh, debug = 1)
      # sink();close(zz);
      # closeAllConnections()
   
      # a1 <- cpp_kNN_rev(a$dist, nh)+1
      # myorder <- order(a$dist[a1], decreasing = FALSE)
      # sort(a1[myorder][1:k])
      
      expect_equal(ist, soll)
   }
})


test_that("bug go to end",{
   skip("")
load("build_ignore/bug_data01.Rda")
   zz <- file("build_ignore/mySink.Rout", open = "wt"); sink(zz); sink(zz, type = "message")
   tmp <- rundtw( Q = Q, C = C, scale = "01", dist_method = "norm1",
                  ws = 20, threshold = 5, lower_bound = T, k = 5 )
   sink();close(zz);
   closeAllConnections()
   tmp
})


test_that("bug counter", {
   skip(" there is no bug!")
   rw <- function(m) cumsum(rnorm(m))
   noise <- function(x) rnorm(length(x))
   deform <- function(x, p) {
      ( simulate_timewarp(x, p, preserve_length = T) + 
           rnorm(1, 0, 10) ) * abs(rnorm(1, 0, 10)) + noise(x)
   }
   infix <- function(x, y, i) {
      if(length(i) == 1){
         new_x <- c(x[1:(i-1)], deform(y, 0.1), x[i:length(x)])
      }else{
         new_x <- infix(infix(x, y, i[1]), y,i[-1])# + length(y))
      }
      return(new_x)
   }
   
   nQ <- 100
   nC <- 1000
   Q <- rw(nQ)
   C0 <- rw(nC)
   ii <- c(30, 550,700)
   C <- infix(C0, Q, ii)
   tmp <- rundtw( Q = Q, C = C, dist_method = "norm1", ws = 10, threshold = NULL,
                          lower_bound = F, k = 0 )
   plot(C, type="l")
   C2 <- C
   C2[-as.vector(sapply(ii, function(i){i:(i+nQ)}))] <- NA
   lines(C2, col="red")
   plot(tmp$dist)

})

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

test_that("figure", {
   skip("")
   dm <- "norm1"
   sp <- "symmetric1"
   WS <- NULL
   noise <- function(nr) cumsum(rnorm(nr))
   
   nx <- 500
   nh <- 100
   nfits <- 2
   
   nn <- nx - nfits * nh# nnoise
   nn <- nn/nfits
   
   h <- noise(nh)
   hscale <- IncDTW::scale( h , type="01")
   x <- noise(0)
   for(i in 1:nfits){
      x <- c(x, noise(nn), h)
   }
   
   a <- rundtw(h, x, k = NULL)
   b <- rundtw(h, x, k = 3)
   
   par(mfrow=c(1,1))
   plot(a$dist)
   points(b$dist, col="red")
   
   
   a <- rundtw(h, x, k = 10)
   c0 <- rundtw(h, x, k = 10, overlap_tol = round(nh/2))
   plot(a$dist)
   points(c0$dist, col="red")

})



test_that("run time univ kNN", {
   skip("runtime comparison")
   
   dm <- "norm1"
   sp <- "symmetric1"
   noise <- function(i) cumsum(rnorm(i))
   # deform <- function(x) (x + rnorm(1, 0, 100)) * abs(rnorm(1, 0, 100))
   
   # speed comparison
   foo_2vec <- function(h, x, ws, k){
      hscale <- IncDTW::scale( h , type="01")
      dis <- sapply(1:(length(x)-length(hscale)+1), function(i){
         y <- IncDTW::scale(x[i:(i+length(hscale)-1)], type = "01")
         IncDTW::dtw2vec(y, hscale, dist_method = dm,
                         step_pattern = sp, ws = ws)$dist
      })
      
      nh <- length(hscale)
      dd <- dis
      best_i <- integer()
      for(i in 1:k){
         ix <- which.min(dd)
         if(dd[ix] < Inf){
            best_i <- c(best_i, ix)
            dd[max(1, (1 + ix - nh)):min((ix + nh - 1), length(dd))] <- Inf   
         }
      }
      return(list(all_values = dis, indices=best_i))
   }
   
   foo_cm <- function(h, x, ws, k){
      hscale <- IncDTW::scale( h , type="01")
      dis <- sapply(1:(length(x)-length(hscale)+1), function(i){
         y <- IncDTW::scale(x[i:(i+length(hscale)-1)], type = "01")
         cm_tmp <- IncDTW::cm(y, hscale, dist_method = dm, ws = ws)
         IncDTW::dtw2vec(cm_tmp, C = "cm", step_pattern = sp, ws = ws)$dist
      })
      
      nh <- length(hscale)
      dd <- dis
      best_i <- integer()
      for(i in 1:k){
         ix <- which.min(dd)
         if(dd[ix] < Inf){
            best_i <- c(best_i, ix)
            dd[max(1, (1 + ix - nh)):min((ix + nh - 1), length(dd))] <- Inf   
         }
      }
      return(list(all_values = dis, indices=best_i))
   }
   
   
   foo_run <- function(h, x, ws, k){
      dis <- rundtw(h, x, dm, sp, k=0, scale = "01", WS, Inf, T)$dist
      dis[is.nan(dis)] <- Inf
      
      nh <- length(hscale)
      dd <- dis
      best_i <- integer()
      for(i in 1:k){
         ix <- which.min(dd)
         if(dd[ix] < Inf){
            best_i <- c(best_i, ix)
            dd[max(1, (1 + ix - nh)):min((ix + nh - 1), length(dd))] <- Inf   
         }
      }
      return(list(all_values = dis, indices=best_i))
   }

         
   nx <- 1000
   nh <- 50
   nfits0 <- 3
   nfits1 <-  min(nfits0, floor(nx/nh))
   
   nfits <- sample(nfits1, size = 1)
   
   nn <- nx - nfits * nh# nnoise
   nn <- nn/nfits
   
   h <- noise(nh)
   x <- c()
   for(i in 1:nfits){
      x <- c(x, noise(nn), h)
   }
   
   WS <- 30
   k <- 5
   mic <- microbenchmark::microbenchmark(foo_cm(h, x, WS,k=k),
                                         foo_2vec(h, x, WS, k=k),
                                         foo_run(h, x, WS, k=k),
                                         rundtw(h, x, dm, sp, k, "01", WS, Inf, T),
                                         times = 30)
   df0 <- maxaR::rel_microbenchmark(mic, 
                                    cols = c("expr", "min", "mean", "median", "max", "neval"), 
                                    rel_expr = "rundtw(h, x, dm, sp, k, '01', WS, Inf, T)")$df
   df0   
   
})


test_that("run time univ", {
   skip("runtime comparison")
   
   dm <- "norm1"
   sp <- "symmetric1"
   nc <- 1
   WS <- NULL
   noise <- function(i) matrix(cumsum(rnorm(i*nc)), nrow =i, ncol=nc)
   # deform <- function(x) (x + rnorm(1, 0, 100)) * abs(rnorm(1, 0, 100))
   
   # speed comparison
   foo_2vec <- function(hscale, x, ws){
      sapply(1:(length(x)-length(hscale)+1), function(i){
         y <- IncDTW::scale(x[i:(i+length(hscale)-1)], type = "01")
         IncDTW::dtw2vec(y, hscale, dist_method = dm,
                         step_pattern = sp, ws = ws)$dist
      })
   }
   
   foo_cm <- function(hscale, x, ws){
      sapply(1:(length(x)-length(hscale)+1), function(i){
         y <- IncDTW::scale(x[i:(i+length(hscale)-1)], type = "01")
         cm_tmp <- IncDTW::cm(y, hscale, dist_method = dm, ws = ws)
         IncDTW::dtw2vec(cm_tmp, C = "cm", step_pattern = sp, ws = ws)$dist
      })
   }
   

   NH <- c(0.05, 0.1, 0.15, 0.2)
   NX <- c(100, 200)
   
   NH <- seq(0.1, 0.5, 0.1)
   NX <- c(100, 200, 300, 500, 1000)
   
   df <- data.frame(expr = character(),
                    min = numeric(), 
                    mean = numeric(), 
                    median = numeric(), 
                    max = numeric(), 
                    neval = integer(), 
                    rel = numeric(), 
                    nh = integer(),
                    nx = integer(),
                    nfits = integer()
                    )
   
   nhh <- NH[1]
   nx <- NX[5]
   
   data.table::fwrite(df, file = "build_ignore/comp_time_test_univ.txt")
   tic <- Sys.time()
   for(nhh in NH){
      print(c(nhh))
      print(Sys.time())
      for(nx in NX){
         print(c(nhh,  nx))
         print(Sys.time())
         
         nh <- nhh*nx
         nfits0 <- 10
         nfits1 <-  min(nfits0, floor(nx/nh))
         
         nfits <- sample(nfits1, size = 1)
         
         nn <- nx - nfits * nh# nnoise
         nn <- nn/nfits
         
         h <- cumsum(rnorm(nh))
         hscale <- IncDTW::scale( h , type="01")
         x <- c()
         for(i in 1:nfits){
            x <- c(x, cumsum(rnorm(nn)), h)
         }
         
         
         
         mic <- microbenchmark::microbenchmark(foo_cm(hscale, x, WS),
                                               foo_2vec(hscale, x, WS),
                                               rundtw(hscale, x, dm, sp, 0, '01', WS, NULL),
                                               rundtw(hscale, x, dm, sp, 0, '01', WS, Inf),
                                               rundtw(hscale, x, dm, sp, 0, "z", WS, NULL),
                                               rundtw(hscale, x, dm, sp, 0, "z", WS, Inf),
                                               rundtw(hscale, x, dm, sp, 0, "none", WS, NULL),
                                               rundtw(hscale, x, dm, sp, 0, "none", WS, Inf),
                                               times = 30)
         df0 <- maxaR::rel_microbenchmark(mic, 
                     cols = c("expr", "min", "mean", "median", "max", "neval"), 
                     rel_expr = "rundtw(hscale, x, dm, sp, 0, \"01\", WS, Inf)")$df
         df0$nh <- nh
         df0$nx <- nx
         df0$nfits <- nfits
         
         data.table::fwrite(df0, file = "build_ignore/comp_time_test_univ.txt", append = T, row.names = F, col.names = F)
         df <- rbind(df, df0)
         
      }
   }
   print(Sys.time())
   
   require(ggplot2)
   
   df <- data.table::fread(file = "build_ignore/comp_time_test_univ.txt")
   df <- df[grep("foo_",df$expr, invert = TRUE), ]
   df <- df[grep("none",df$expr, invert = TRUE), ]
   ggplot(df) + geom_line(aes(x = nh/nx, y = rel, group = expr, col = expr), size = 2)+
      facet_grid(~nx)
   #
   
   
   
   x <- matrix(cumsum(rnorm(3 * 1000)), ncol = 3)
   hscale <- matrix(cumsum(rnorm(3 * 100)), ncol = 3)
   
   
   mic <- microbenchmark::microbenchmark(rundtw(hscale, x, step_pattern = "symmetric1", scale = '01'  , dist_method = "norm2_square", k = 0, ws = NULL, threshold = NULL, lower_bound = FALSE ),
                                         rundtw(hscale, x, step_pattern = "symmetric1", scale = 'z'   , dist_method = "norm2_square", k = 0, ws = NULL, threshold = NULL, lower_bound = FALSE ),
                                         rundtw(hscale, x, step_pattern = "symmetric1", scale = 'none', dist_method = "norm2_square", k = 0, ws = NULL, threshold = NULL, lower_bound = FALSE ),
                                         times = 5)
   levels(mic$expr) <- c("f(01)", "f(z)", "f(none)")
   mic

})







test_that("knn 01-scale vs. knn z-scale", {
   skip("comparison 01-z")
   
   dm <- "norm1"
   sp <- "symmetric1"
   WS <- 10
   
   noise <- function(i) cumsum(rnorm(i))

   h <- c(20, noise(10))
   hscale01 <- IncDTW::scale(h, type="01")
   hscalez <- IncDTW::scale(h, type="z")
   
   par(mfrow=c(3,1))
   plot(h, type="l")
   plot(hscale01, type="l")
   plot(hscalez, type="l")
   par(mfrow=c(1,1))
   
   x <- c(noise(10), h, noise(10), h, noise(10), h, noise(10), h, noise(10))
   
   ret01 <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, 
                 ws = WS, scale = "01", threshold = NULL, k = 10)
   retz  <- rundtw(Q = h, C = x, dist_method = dm, step_pattern = sp, 
                   ws = WS, scale = "z", threshold = NULL, k = 10)
      
   ret01$knn_indices
   retz$knn_indices
   
})






