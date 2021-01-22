context("scaling time series")


test_that("test scale with zero variance", {
   #skip("blabla")
   x <- rep(rnorm(1), 100)
   
   fz  <- function(x){  (x-mean(x))  }
   f01 <- function(x){  (x-min(x)) }
   
   
   ist <- scale(x, type = "z")
   soll <- fz(x)
   expect_equal(ist, soll)
   
   ist <- scale(x, type = "01")
   soll <- f01(x)
   expect_equal(ist, soll)

})


test_that("test scale with low variance", {
   #skip("blabla")
   x <- rep(rnorm(1)/10^6, 100)
   x[1] <- x[1] * 2
   
   fz  <- function(x){  (x-mean(x))/sd(x)  }
   f01 <- function(x){  (x-min(x))/(max(x)-min(x)) }
   
   
   ist <- scale(x, type = "z", threshold = 1e-12)
   soll <- fz(x)
   expect_equal(ist, soll)
   
   ist <- scale(x, type = "01", threshold = 1e-12)
   soll <- f01(x)
   expect_equal(ist, soll)
  
})


test_that("test zscale", {
   #skip("blabla")
   eps <- 1/10^7
   x <- rnorm(100, sample(100,1), sample(10,1))
   
   f0 <- function(x){  (x-mean(x))/sd(x)  }
   f1 <- function(x, mu0, sd0){  (x-mu0)/sd0  }
   
   mm <- rnorm(1)
   ss <- abs(rnorm(1))
   expect_equal(max(abs(f0(x) - scale(x, type = "z"))) < eps, TRUE)
   expect_equal(max(abs(f1(x, mu0 = mm, sd0 = ss) - 
                        scale(x, type = "z", xmean = mm, xsd = ss))) < eps, TRUE)
   expect_equal(max(abs(f0(x) - 
                        scale(x, type = "z", xmean = mean(x)))) < eps, TRUE)
   expect_equal(max(abs(f0(x) - 
                        scale(x, type = "z", xsd = sd(x)))) < eps, TRUE)
   
   
})



test_that("test znorm speed", {
   skip("speed comparison")
   
   x <- rnorm(10^6, sample(100,1), sample(10,1))
   
   f0 <- function(x){  (x-mean(x))/sd(x)  }
   f1 <- function(x, mu0, sd0){  (x-mu0)/sd0  }
   
   mm <- rnorm(1)
   ss <- rnorm(1)
   microbenchmark::microbenchmark(f0(x),
                                  cpp_znorm(x),
                                  cpp_znorm(x, mm, ss))
   
})



test_that("test scale01", {
   #skip("blabla")
   eps <- 1/10^7
   x <- rnorm(100, sample(100,1), sample(10,1))
   f0 <- function(x){  (x-min(x))/(max(x)-min(x))  }
   f1 <- function(x, mi, ma){  (x-mi)/(ma - mi)  }
   
   mmi <- rnorm(1)
   mma <- mmi + abs(rnorm(1))
   
   expect_equal(max(abs(f0(x) - scale(x, type = "01"))) < eps, TRUE)
   expect_equal(max(abs(f1(x, mi = mmi, ma = mma) - 
                           scale(x, type = "01", xmax = mma, xmin = mmi))) < eps, TRUE)
   expect_equal(max(abs(f0(x) - 
                        scale(x, type = "01", xmax = max(x)))) < eps, TRUE)
   expect_equal(max(abs(f0(x) - 
                           scale(x, type = "01", xmin = min(x)))) < eps, TRUE)
   
})



test_that("test scale01 speed", {
   skip("speed comparison")
   
   x <- rnorm(10^4, sample(100,1), sample(10,1))
   
   f0 <- function(x){  (x-min(x))/(max(x)-min(x))  }
   f1 <- function(x, mi, ma){  (x-mi)/(ma - mi)  }
   
   mmi <- rnorm(1)
   mma <- rnorm(1)
   
   microbenchmark::microbenchmark(f0(x),
                                  cpp_norm01(x),
                                  cpp_norm01(x, min_in = mmi, max_in = mma))
   
})
