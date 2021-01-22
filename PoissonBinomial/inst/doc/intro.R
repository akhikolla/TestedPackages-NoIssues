## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo = FALSE------------------------------------------------------
library(PoissonBinomial)

## ----ex-opdb------------------------------------------------------------------
# Case 1
dpbinom(NULL, rep(0.3, 7))
dbinom(0:7, 7, 0.3) # equal results

dpbinom(NULL, c(0, 0, 0, 0, 0, 0, 0)) # only 0 is observable
dpbinom(0, c(0, 0, 0, 0, 0, 0, 0)) # confirmation

dpbinom(NULL, c(1, 1, 1, 1, 1, 1, 1)) # only 7 is observable
dpbinom(7, c(1, 1, 1, 1, 1, 1, 1)) # confirmation

# Case 2
dpbinom(NULL, c(0, 0, 0, 0, 1, 1, 1)) # only 3 is observable
dpbinom(3, c(0, 0, 0, 0, 1, 1, 1)) # confirmation

# Case 3
dpbinom(NULL, c(0, 0, 0.1, 0.2, 0.4, 0.8, 1)) # only 1-5 are observable
dpbinom(1:5, c(0, 0, 0.1, 0.2, 0.4, 0.8, 1)) # confirmation

dpbinom(NULL, c(0, 0, 0.4, 1)) # only 1 and 2 are observable
dpbinom(1:2, c(0, 0, 0.4, 1)) # confirmation

## ----ex-gpdb------------------------------------------------------------------
set.seed(1)
pp <- runif(7)
va <- sample(0:6, 7, TRUE)
vb <- sample(0:6, 7, TRUE)

# Case 1
dgpbinom(NULL, pp, rep(1, 7), rep(0, 7))
dpbinom(NULL, pp) # equal results

dgpbinom(NULL, pp, rep(0, 7), rep(1, 7))
dpbinom(NULL, 1 - pp) # equal results

dgpbinom(NULL, pp, c(rep(1, 3), rep(0, 4)), c(rep(0, 3), rep(1, 4)))
dpbinom(NULL, c(pp[1:3], 1 - pp[4:7])) # reorder for 0 and 1; equal results

# Case 2 a)
dgpbinom(NULL, rep(0, 7), rep(4, 7), rep(2, 7)) # only 14 is observable
dgpbinom(7 * 2, rep(0, 7), rep(4, 7), rep(2, 7)) # confirmation

# Case 2 b)
dgpbinom(NULL, rep(1, 7), rep(4, 7), rep(2, 7)) # only 28 is observable
dgpbinom(7 * 4, rep(1, 7), rep(4, 7), rep(2, 7)) # confirmation

# Case 2 c)
dgpbinom(NULL, rep(0.3, 7), rep(4, 7), rep(2, 7))
dbinom(0:7, 7, 0.3) # equal results, but on different support set

# Case 2 d)
dgpbinom(NULL, pp, rep(4, 7), rep(2, 7))
dpbinom(NULL, pp) # equal results, but on different support set

# Case 3 a)
dgpbinom(NULL, c(0, 0, 0, 0, 0, 0, 0), va, vb) # only sum(vb) is observable
dgpbinom(sum(vb), rep(0, 7), va, vb) # confirmation

# Case 3 b)
dgpbinom(NULL, c(1, 1, 1, 1, 1, 1, 1), va, vb) # only sum(va) is observable
dgpbinom(sum(va), rep(1, 7), va, vb) # confirmation

# Case 3 c)
dgpbinom(NULL, c(0, 0, 0, 1, 1, 1, 1), va, vb) # only sum(va[4:7], vb[1:3]) is observable
dgpbinom(sum(va[4:7], vb[1:3]), c(0, 0, 0, 1, 1, 1, 1), va, vb) # confirmation

# Case 4
dgpbinom(NULL, c(0, 0, 0.3, 0.6, 1, 1, 1), va, vb)
sure <- sum(va[5:7], vb[1:2])
x.transf <- sum(pmin(va[3:4], vb[3:4])):sum(pmax(va[3:4], vb[3:4]))
dgpbinom(sure + x.transf, c(0, 0, 0.3, 0.6, 1, 1, 1), va, vb)
dgpbinom(x.transf, c(0.3, 0.6), va[3:4], vb[3:4]) # equal results

dgpbinom(NULL, c(0, 0, 0, 0.6, 1, 1, 1), va, vb)
sure <- sum(va[5:7], vb[1:3])
x.transf <- va[4]:vb[4]
dgpbinom(sure + x.transf, c(0, 0, 0, 0.6, 1, 1, 1), va, vb)
dgpbinom(x.transf, 0.6, va[4], vb[4]) # equal results; essentially transformed Bernoulli

