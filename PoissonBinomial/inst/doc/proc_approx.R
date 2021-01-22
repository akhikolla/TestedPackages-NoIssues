## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo = FALSE------------------------------------------------------
library(PoissonBinomial)

## ----pa1----------------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)

dpbinom(NULL, pp, wt, "Poisson")
ppbinom(NULL, pp, wt, "Poisson")

## ----pa2----------------------------------------------------------------------
set.seed(1)

# U(0, 1) random probabilities of success
pp <- runif(20)
dpbinom(NULL, pp, method = "Poisson")
dpbinom(NULL, pp)
summary(dpbinom(NULL, pp, method = "Poisson") - dpbinom(NULL, pp))

# U(0, 0.01) random probabilities of success
pp <- runif(20, 0, 0.01)
dpbinom(NULL, pp, method = "Poisson")
dpbinom(NULL, pp)
summary(dpbinom(NULL, pp, method = "Poisson") - dpbinom(NULL, pp))

## ----am1----------------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)
mean(rep(pp, wt))

dpbinom(NULL, pp, wt, "Mean")
ppbinom(NULL, pp, wt, "Mean")

## ----am2----------------------------------------------------------------------
set.seed(1)

# U(0, 1) random probabilities of success
pp <- runif(20)
dpbinom(NULL, pp, method = "Mean")
dpbinom(NULL, pp)
summary(dpbinom(NULL, pp, method = "Mean") - dpbinom(NULL, pp))

# U(0.3, 0.5) random probabilities of success
pp <- runif(20, 0.3, 0.5)
dpbinom(NULL, pp, method = "Mean")
dpbinom(NULL, pp)
summary(dpbinom(NULL, pp, method = "Mean") - dpbinom(NULL, pp))

# U(0.39, 0.41) random probabilities of success
pp <- runif(20, 0.39, 0.41)
dpbinom(NULL, pp, method = "Mean")
dpbinom(NULL, pp)
summary(dpbinom(NULL, pp, method = "Mean") - dpbinom(NULL, pp))

## ----gma1---------------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)
prod(rep(pp, wt))^(1/sum(wt))

dpbinom(NULL, pp, wt, "GeoMean")
ppbinom(NULL, pp, wt, "GeoMean")

## ----gma2---------------------------------------------------------------------
set.seed(1)

# U(0, 1) random probabilities of success
pp <- runif(20)
dpbinom(NULL, pp, method = "GeoMean")
dpbinom(NULL, pp)
summary(dpbinom(NULL, pp, method = "GeoMean") - dpbinom(NULL, pp))

# U(0.4, 0.6) random probabilities of success
pp <- runif(20, 0.4, 0.6)
dpbinom(NULL, pp, method = "GeoMean")
dpbinom(NULL, pp)
summary(dpbinom(NULL, pp, method = "GeoMean") - dpbinom(NULL, pp))

# U(0.49, 0.51) random probabilities of success
pp <- runif(20, 0.49, 0.51)
dpbinom(NULL, pp, method = "GeoMean")
dpbinom(NULL, pp)
summary(dpbinom(NULL, pp, method = "GeoMean") - dpbinom(NULL, pp))

## ----gmb1---------------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)
1 - prod(1 - rep(pp, wt))^(1/sum(wt))

dpbinom(NULL, pp, wt, "GeoMeanCounter")
ppbinom(NULL, pp, wt, "GeoMeanCounter")

## ----gmb2---------------------------------------------------------------------
set.seed(1)

# U(0, 1) random probabilities of success
pp <- runif(20)
dpbinom(NULL, pp, method = "GeoMeanCounter")
dpbinom(NULL, pp)
summary(dpbinom(NULL, pp, method = "GeoMeanCounter") - dpbinom(NULL, pp))

# U(0.4, 0.6) random probabilities of success
pp <- runif(20, 0.4, 0.6)
dpbinom(NULL, pp, method = "GeoMeanCounter")
dpbinom(NULL, pp)
summary(dpbinom(NULL, pp, method = "GeoMeanCounter") - dpbinom(NULL, pp))

# U(0.49, 0.51) random probabilities of success
pp <- runif(20, 0.49, 0.51)
dpbinom(NULL, pp, method = "GeoMeanCounter")
dpbinom(NULL, pp)
summary(dpbinom(NULL, pp, method = "GeoMeanCounter") - dpbinom(NULL, pp))

## ----na1-ord------------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)

dpbinom(NULL, pp, wt, "Normal")
ppbinom(NULL, pp, wt, "Normal")

## ----na2-ord------------------------------------------------------------------
set.seed(1)

# 10 random probabilities of success
pp <- runif(10)
dpn <- dpbinom(NULL, pp, method = "Normal")
dpd <- dpbinom(NULL, pp)
idx <- which(dpn != 0 & dpd != 0)
summary((dpn - dpd)[idx])

# 1000 random probabilities of success
pp <- runif(1000)
dpn <- dpbinom(NULL, pp, method = "Normal")
dpd <- dpbinom(NULL, pp)
idx <- which(dpn != 0 & dpd != 0)
summary((dpn - dpd)[idx])

# 100000 random probabilities of success
pp <- runif(100000)
dpn <- dpbinom(NULL, pp, method = "Normal")
dpd <- dpbinom(NULL, pp)
idx <- which(dpn != 0 & dpd != 0)
summary((dpn - dpd)[idx])

## ----rna1-ord-----------------------------------------------------------------
set.seed(1)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)

dpbinom(NULL, pp, wt, "RefinedNormal")
ppbinom(NULL, pp, wt, "RefinedNormal")

## ----rna2-ord-----------------------------------------------------------------
set.seed(1)

# 10 random probabilities of success
pp <- runif(10)
dpn <- dpbinom(NULL, pp, method = "RefinedNormal")
dpd <- dpbinom(NULL, pp)
idx <- which(dpn != 0 & dpd != 0)
summary((dpn - dpd)[idx])

# 1000 random probabilities of success
pp <- runif(1000)
dpn <- dpbinom(NULL, pp, method = "RefinedNormal")
dpd <- dpbinom(NULL, pp)
idx <- which(dpn != 0 & dpd != 0)
summary((dpn - dpd)[idx])

# 100000 random probabilities of success
pp <- runif(100000)
dpn <- dpbinom(NULL, pp, method = "RefinedNormal")
dpd <- dpbinom(NULL, pp)
idx <- which(dpn != 0 & dpd != 0)
summary((dpn - dpd)[idx])

## ----benchmark-ord------------------------------------------------------------
library(microbenchmark)
set.seed(1)

f1 <- function() dpbinom(NULL, runif(4000), method = "Normal")
f2 <- function() dpbinom(NULL, runif(4000), method = "RefinedNormal")
f3 <- function() dpbinom(NULL, runif(4000), method = "Poisson")
f4 <- function() dpbinom(NULL, runif(4000), method = "Mean")
f5 <- function() dpbinom(NULL, runif(4000), method = "GeoMean")
f6 <- function() dpbinom(NULL, runif(4000), method = "GeoMeanCounter")
f7 <- function() dpbinom(NULL, runif(4000), method = "DivideFFT")

microbenchmark(f1(), f2(), f3(), f4(), f5(), f6(), f7())

## ----na1-gen------------------------------------------------------------------
set.seed(2)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)
va <- sample(0:10, 10, TRUE)
vb <- sample(0:10, 10, TRUE)

dgpbinom(NULL, pp, va, vb, wt, "Normal")
pgpbinom(NULL, pp, va, vb, wt, "Normal")

## ----na2-gen------------------------------------------------------------------
set.seed(2)

# 10 random probabilities of success
pp <- runif(10)
va <- sample(0:10, 10, TRUE)
vb <- sample(0:10, 10, TRUE)
dpn <- dgpbinom(NULL, pp, va, vb, method = "Normal")
dpd <- dgpbinom(NULL, pp, va, vb)
idx <- which(dpn != 0 & dpd != 0)
summary((dpn - dpd)[idx])

# 100 random probabilities of success
pp <- runif(100)
va <- sample(0:100, 100, TRUE)
vb <- sample(0:100, 100, TRUE)
dpn <- dgpbinom(NULL, pp, va, vb, method = "Normal")
dpd <- dgpbinom(NULL, pp, va, vb)
idx <- which(dpn != 0 & dpd != 0)
summary((dpn - dpd)[idx])

# 1000 random probabilities of success
pp <- runif(1000)
va <- sample(0:1000, 1000, TRUE)
vb <- sample(0:1000, 1000, TRUE)
dpn <- dgpbinom(NULL, pp, va, vb, method = "Normal")
dpd <- dgpbinom(NULL, pp, va, vb)
idx <- which(dpn != 0 & dpd != 0)
summary((dpn - dpd)[idx])

## ----rna1-gen-----------------------------------------------------------------
set.seed(2)
pp <- runif(10)
wt <- sample(1:10, 10, TRUE)
va <- sample(0:10, 10, TRUE)
vb <- sample(0:10, 10, TRUE)
dgpbinom(NULL, pp, va, vb, wt, "RefinedNormal")
pgpbinom(NULL, pp, va, vb, wt, "RefinedNormal")

## ----rna2-gen-----------------------------------------------------------------
set.seed(2)

# 10 random probabilities of success
pp <- runif(10)
va <- sample(0:10, 10, TRUE)
vb <- sample(0:10, 10, TRUE)
dpn <- dgpbinom(NULL, pp, va, vb, method = "RefinedNormal")
dpd <- dgpbinom(NULL, pp, va, vb)
idx <- which(dpn != 0 & dpd != 0)
summary((dpn - dpd)[idx])

# 100 random probabilities of success
pp <- runif(100)
va <- sample(0:100, 100, TRUE)
vb <- sample(0:100, 100, TRUE)
dpn <- dgpbinom(NULL, pp, va, vb, method = "RefinedNormal")
dpd <- dgpbinom(NULL, pp, va, vb)
idx <- which(dpn != 0 & dpd != 0)
summary((dpn - dpd)[idx])

# 1000 random probabilities of success
pp <- runif(1000)
va <- sample(0:1000, 1000, TRUE)
vb <- sample(0:1000, 1000, TRUE)
dpn <- dgpbinom(NULL, pp, va, vb, method = "RefinedNormal")
dpd <- dgpbinom(NULL, pp, va, vb)
idx <- which(dpn != 0 & dpd != 0)
summary((dpn - dpd)[idx])

## ----benchmark-gen------------------------------------------------------------
library(microbenchmark)
n <- 1500
set.seed(2)
va <- sample(1:50, n, TRUE)
vb <- sample(1:50, n, TRUE)

f1 <- function() dgpbinom(NULL, runif(n), va, vb, method = "Normal")
f2 <- function() dgpbinom(NULL, runif(n), va, vb, method = "RefinedNormal")
f3 <- function() dgpbinom(NULL, runif(n), va, vb, method = "DivideFFT")

microbenchmark(f1(), f2(), f3())

