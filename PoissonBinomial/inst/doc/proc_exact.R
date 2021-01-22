## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo = FALSE------------------------------------------------------
library(PoissonBinomial)

## ----directconv-ord-----------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)

dpbinom(NULL, pp, wt, "Convolve")
ppbinom(NULL, pp, wt, "Convolve")

## ----divide1-ord--------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)

dpbinom(NULL, pp, wt, "DivideFFT")
ppbinom(NULL, pp, wt, "DivideFFT")

## ----divide2-ord--------------------------------------------------------------
set.seed(1)
pp1 <- runif(751)
pp2 <- pp1[1:750]

sum(abs(dpbinom(NULL, pp2, method = "DivideFFT") - dpbinom(NULL, pp2, method = "Convolve")))
sum(abs(dpbinom(NULL, pp1, method = "DivideFFT") - dpbinom(NULL, pp1, method = "Convolve")))

## ----divide3-ord--------------------------------------------------------------
set.seed(1)
pp1 <- runif(751)

d1 <- dpbinom(NULL, pp1, method = "DivideFFT")
d2 <- dpbinom(NULL, pp1, method = "Convolve")

min(d1[d1 > 0])
min(d2[d2 > 0])

## ----dftcf-ord----------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)

dpbinom(NULL, pp, wt, "Characteristic")
ppbinom(NULL, pp, wt, "Characteristic")

## ----rf1-ord------------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)

dpbinom(NULL, pp, wt, "Recursive")
ppbinom(NULL, pp, wt, "Recursive")

## ----rf2-ord------------------------------------------------------------------
set.seed(1)
pp <- runif(1000)
wt <- sample(1:10, 1000, TRUE)

sum(abs(dpbinom(NULL, pp, wt, "Convolve") - dpbinom(NULL, pp, wt, "Recursive")))

## ----benchmark-ord------------------------------------------------------------
library(microbenchmark)
set.seed(1)

f1 <- function() dpbinom(NULL, runif(5000), method = "DivideFFT")
f2 <- function() dpbinom(NULL, runif(5000), method = "Convolve")
f3 <- function() dpbinom(NULL, runif(5000), method = "Recursive")
f4 <- function() dpbinom(NULL, runif(5000), method = "Characteristic")

microbenchmark(f1(), f2(), f3(), f4())

## ----directconv-gen-----------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)
va <- sample(0:10, 10, TRUE)
vb <- sample(0:10, 10, TRUE)

dgpbinom(NULL, pp, va, vb, wt, "Convolve")
pgpbinom(NULL, pp, va, vb, wt, "Convolve")

## ----divide1-gen--------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)
va <- sample(0:10, 10, TRUE)
vb <- sample(0:10, 10, TRUE)

dgpbinom(NULL, pp, va, vb, wt, "DivideFFT")
pgpbinom(NULL, pp, va, vb, wt, "DivideFFT")

## ----divide2-gen--------------------------------------------------------------
set.seed(1)
pp1 <- runif(250)
va1 <- sample(0:50, 250, TRUE)
vb1 <- sample(0:50, 250, TRUE)
pp2 <- pp1[1:248]
va2 <- va1[1:248]
vb2 <- vb1[1:248]

sum(abs(dgpbinom(NULL, pp1, va1, vb1, method = "DivideFFT")
        - dgpbinom(NULL, pp1, va1, vb1, method = "Convolve")))

sum(abs(dgpbinom(NULL, pp2, va2, vb2, method = "DivideFFT")
        - dgpbinom(NULL, pp2, va2, vb2, method = "Convolve")))

## ----divide3-gen--------------------------------------------------------------
d1 <- dgpbinom(NULL, pp1, va1, vb1, method = "DivideFFT")
d2 <- dgpbinom(NULL, pp1, va1, vb1, method = "Convolve")

min(d1[d1 > 0])
min(d2[d2 > 0])

## ----dftcf-gen----------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)
va <- sample(0:10, 10, TRUE)
vb <- sample(0:10, 10, TRUE)

dgpbinom(NULL, pp, va, vb, wt, "Characteristic")
pgpbinom(NULL, pp, va, vb, wt, "Characteristic")

## ----benchmark-gen------------------------------------------------------------
library(microbenchmark)
n <- 2500
set.seed(1)
va <- sample(1:50, n, TRUE)
vb <- sample(1:50, n, TRUE)

f1 <- function() dgpbinom(NULL, runif(n), va, vb, method = "DivideFFT")
f2 <- function() dgpbinom(NULL, runif(n), va, vb, method = "Convolve")
f3 <- function() dgpbinom(NULL, runif(n), va, vb, method = "Characteristic")

microbenchmark(f1(), f2(), f3())

