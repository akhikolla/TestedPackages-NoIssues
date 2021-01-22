context("test find_peaks")


test_that("min equal expected", {
   
   # finds also the border peak 29
   x = c(1:10, 9:1, 2:11)
   ist <-  find_peaks(x, w=3, get_min=TRUE)
   soll <- c(1, 19)
   expect_equal(ist, soll)
   
   # what means the neigbohood parameter m:
   # at least m-many values need to be inbetween two consecutive peaks
   x = c(1:10, 9, 9, 11, 9:8, 7)
   ist <-  find_peaks(x, w=3, get_min=TRUE)
   soll <- c(1, 16)
   expect_equal(ist, soll)
   
   x = c(1:10, 9, 9,9, 11, 9:8, 7)
   ist <-  find_peaks(x, w=3, get_min=TRUE, strict = TRUE)
   soll <- c(1,  17)
   expect_equal(ist, soll)
   
   x = c(1:10, 9, 9,9, 11, 9:8, 7)
   ist <-  find_peaks(x, w=3, get_min=TRUE, strict = FALSE)
   soll <- c(1, 12, 17)
   expect_equal(ist, soll)
   
})


test_that("saddle points equal expected", {
   # also saddle points
   x = c(1:10, rep(10, 10), 10:1)
   ist <-  find_peaks(x, w=3, get_min=FALSE, strict = TRUE)
   soll <- integer()
   expect_equal(ist, soll)

   x = c(1:10, rep(10, 10), 10:1)
   ist <-  find_peaks(x, w=3, get_min=FALSE, strict = FALSE)
   soll <- 10:21
   expect_equal(ist, soll)
   
   
   x = c(1:10, rep(10, 10), 10:21)
   ist <-  find_peaks(x, w=3, get_min=FALSE, strict = FALSE)
   soll <- c(10:18, 32)
   expect_equal(ist, soll)
   
})


test_that("nan values", {
   # also saddle points
   
   x <-  c(1 ,rep(c(NaN), 10))
   ist <-  find_peaks(x, w=3, get_min=TRUE, strict = FALSE)
   soll <- 1
   expect_equal(ist, soll)
   
   x <-  c(1 ,rep(c(NA), 10))
   ist <-  find_peaks(x, w=3, get_min=TRUE, strict = FALSE)
   soll <- 1
   expect_equal(ist, soll)
   
   x <-  c(1 ,rep(c(NaN), 10))
   ist <-  find_peaks(x, w=3, get_min=FALSE, strict = FALSE)
   soll <- 1
   expect_equal(ist, soll)
   
   x <-  c(1 ,rep(c(NA), 10))
   ist <-  find_peaks(x, w=3, get_min=FALSE, strict = FALSE)
   soll <- 1
   expect_equal(ist, soll)
   
})

test_that("max equal expected", {
  
   # finds also the border peak 29
   x = c(1:10, 9:1, 2:11)
   ist <-  find_peaks(x, w=3, get_min=FALSE)
   soll <- c(10, 29)
   expect_equal(ist, soll)
   
   # what means the neigbohood parameter m:
   # at least m-many values need to be inbetween two consecutive peaks
   x = c(1:10, 9, 9, 11, 9:8, 7)
   ist <-  find_peaks(x, w=3, get_min=FALSE)
   soll <- c(13)
   expect_equal(ist, soll)
   
   x = c(1:10, 9, 9,9, 11, 9:8, 7)
   ist <-  find_peaks(x, w=3, get_min=FALSE)
   soll <- c(10, 14)
   expect_equal(ist, soll)
   
   
})

