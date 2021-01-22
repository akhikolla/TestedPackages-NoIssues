## ---- include=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)
library("distr6")
set.seed(42)

## -----------------------------------------------------------------------------
Normal$new()

## -----------------------------------------------------------------------------
Normal$new(mean = 2, sd = 2)
Normal$new(mean = 3, prec = 0.5)

## -----------------------------------------------------------------------------
N <- Normal$new()
N$print()
N$parameters()

## -----------------------------------------------------------------------------
N$setParameterValue(prec = 2)
N$getParameterValue("prec")

## -----------------------------------------------------------------------------
N$parameters()

## -----------------------------------------------------------------------------
N$parameters()$print(hide_cols = NULL)

## -----------------------------------------------------------------------------
N$setParameterValue(var = 3)$getParameterValue("var")

## -----------------------------------------------------------------------------
N$print()

## -----------------------------------------------------------------------------
N$summary()
N$summary(full = F)

## -----------------------------------------------------------------------------
N$properties
N$traits

## -----------------------------------------------------------------------------
N$pdf(2) # dnorm(2)
N$cdf(2) # pnorm(2)
N$quantile(0.42) # qnorm(2)
N$rand(2) # rnorm(2)

## -----------------------------------------------------------------------------
B <- Beta$new(shape1 = 0.582, shape2 = 1.2490)
B$pdf(2) # dbeta(2, 0.582, 1.2490)
B$cdf(2) # pbeta(2, 0.582, 1.2490)
B$quantile(0.42) # qbeta(2, 0.582, 1.2490)
B$rand(2) # rbeta(2, 0.582, 1.2490)

## -----------------------------------------------------------------------------
N$cdf(3, lower.tail = FALSE, log.p = TRUE) == pnorm(3, lower.tail = FALSE, log.p = TRUE)

## -----------------------------------------------------------------------------
N$mean()
N$variance()
N$entropy() # Note default is base 2
N$mgf(2)
N$cf(1)

## ---- eval=FALSE, echo=FALSE, message=FALSE, warning=FALSE--------------------
#  ?Normal

## -----------------------------------------------------------------------------
head(listDistributions())
head(listDistributions(simplify = TRUE))
# Lists discrete distributions only
head(listDistributions(filter = list(valuesupport = "discrete")))

# Multiple filters can be used, note this is case-insensitive
head(listDistributions(filter = list(VaLueSupport = "continuous", package = "extraDistr")))

## -----------------------------------------------------------------------------
library(magrittr)
N$print()
print(N)
N %>% print()

N$pdf(2)
pdf(N, 2)
N %>% pdf(2)

