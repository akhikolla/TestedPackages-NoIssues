## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
  )

## ------------------------------------------------------------------------
### uncomment to install packages:
# install.packages("devtools", "coda", "RcppArmadillo", "RcppProgress", "Rglpk", "quadprog")
# devtools::install_github("danheck/multinomineq")

library("multinomineq")
library("coda")
set.seed(1234)

## ------------------------------------------------------------------------
# total number of observations for each paired comparison:
n <- 25  
# Frequencies of choosing "Option H" in Description condition:
HER_desc <- c(9, 16, 16, 7, 12, 16)
# Frequencies of choosing "Option H" in Experience condition:
HER_exp  <- c(22, 11, 7, 14, 5, 3)

## ------------------------------------------------------------------------
# Underweighting of small probabilities CPT:
# (1 = choose "Option H"; 0 = do not choose "Option H")
V <- matrix(c(0, 0, 0, 0, 0, 0,  # first predicted pattern
              0, 0, 0, 1, 0, 0,  # second predicted pattern
              0, 0, 1, 1, 0, 0,  # ...
              1, 0, 0, 0, 0, 0,
              1, 0, 0, 1, 0, 0,
              1, 0, 1, 1, 0, 0,
              1, 1, 0, 0, 0, 0,
              1, 1, 0, 0, 1, 0,
              1, 1, 0, 1, 1, 0,
              1, 1, 0, 0, 1, 1,
              1, 1, 0, 1, 0, 0,
              1, 1, 0, 1, 1, 1,
              1, 1, 1, 1, 0, 0,
              1, 1, 1, 1, 1, 0,
              1, 1, 1, 1, 1, 1), 
            ncol = 6, byrow = TRUE)

## ------------------------------------------------------------------------
# define a specific vector of choice probabilities:
p_observed <- c(.25, .48, .93, .10, .32, .50)

# check whether p is inside the convex hull defined by V:
inside(p_observed, V = V)

# to find a (random) probability that is inside the convex hull:
find_inside(V = V, random = TRUE)

## ------------------------------------------------------------------------
A <- matrix(c( 0,   0,   -1,    0,    0,    0,  # first inequality
               0,   0,    0,    0,    0,   -1,  # second inequality
              -1,   1,    0,    0,    0,    0,  # ...
               0,  -1,    0,    0,    1,    0,
               0,   0,    0,    0,   -1,    1,
               0,   0,    1,   -1,    0,    0,
               0,   0,    0,    1,    0,    0,
               1,   0,    0,    0,    0,    0), 
            ncol = 6, byrow = TRUE)

b <- c(0, 0, 0, 0, 0, 0, 1, 1)

## ------------------------------------------------------------------------
p_observed <- c(.25, .48, .93, .10, .32, .50)
all(A %*% p_observed <= b)

# corresponding function in multinomineq:
inside(p_observed, A, b)

# find a point that satisfies the constraints:
find_inside(A, b, random = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  ### (not run):
#  devtools::install_github("TasCL/rPorta")
#  Ab <- V_to_Ab(V)
#  Ab

## ---- eval = FALSE-------------------------------------------------------
#  # fit data from Description and Experience condition:
#  p_D <- sampling_binom(k = HER_desc, n = n, V = V, M = 2000, progress = FALSE)
#  p_E <- sampling_binom(k = HER_exp, n = n, V = V, M = 2000, progress = FALSE)

## ------------------------------------------------------------------------
# fit data from Description and Experience condition:
p_D <- sampling_binom(k = HER_desc, n = n, A = A, b = b, 
                      M = 2000, progress = FALSE)
p_E <- sampling_binom(k = HER_exp, n = n, A = A, b = b, 
                      M = 2000, progress = FALSE)

## ---- fig.width=7, fig.height=5.5----------------------------------------
# check convergence of MCMC sampler (should look like hairy caterpillars):
plot(p_D)
# summarize posterior samples:
summary(p_D)

## ------------------------------------------------------------------------
# posterior-predictive p-value for Description condition:
ppp_binom(prob = p_D, k = HER_desc, n = n)
ppp_binom(prob = p_D, k = HER_exp, n = n)

## ----eval = FALSE--------------------------------------------------------
#  # compute Bayes factor using the V-representation
#  bf_D <- bf_binom(k = HER_desc, n = n, V = V, M = 10000, progress = FALSE)
#  bf_E <- bf_binom(k = HER_exp,  n = n, V = V, M = 10000, progress = FALSE)

## ------------------------------------------------------------------------
# compute Bayes factor using the Ab-representation
bf_D <- bf_binom(HER_desc, n = n, A = A, b = b, M = 10000, progress = FALSE)
bf_E <- bf_binom(HER_exp,  n = n, A = A, b = b, M = 10000, progress = FALSE)

bf_D
bf_E

## ------------------------------------------------------------------------
# proportion of *prior* samples satisfying inequality constraints
cnt <- count_binom(k = 0, n = 0, A = A, b = b, M = 50000, progress = FALSE)
cnt

# proportion of *prior* samples satisfying inequality constraints 
# (dataset: Experience)
cnt_E <- count_binom(k = HER_exp, n = n, A = A, b = b, M = 10000, progress = FALSE)
cnt_E

# obtain Bayes factor (including error of approximation)
count_to_bf(posterior = cnt_E, prior = cnt)

## ------------------------------------------------------------------------
# use a stepwise counting approach for the "Description" condition:
# ("cmin" = minimum of samples that satisfy constraints before sampling stops)
cnt_D <- count_binom(k = HER_desc, n = n, A = A, b = b, 
                     M = 5000, cmin = 1, steps = c(3, 5, 7), progress = FALSE)
cnt_D

# obtain Bayes factor (including error of approximation)
count_to_bf(posterior = cnt_D, prior = cnt)

## ------------------------------------------------------------------------
# binomial data format
n <- 25  
HER_desc <- c(9, 16, 16, 7, 12, 16)

# multinomial data format  ("k" must be a vector!)
k_multinom <- c(9, 16,   # first binary gamble
                16, 9,   # second binary gamble
                16, 9,   # ...
                7, 17,   
                12, 13,   
                16, 9)
options <- c(2, 2, 2, 2, 2, 2)  # 2 options for each type

## ------------------------------------------------------------------------
mn <- binom_to_multinom(HER_desc, n)
mn

## ------------------------------------------------------------------------
posterior <- sampling_multinom(k = mn$k, options = mn$options, A = A, b = b)
bayesfactor <- bf_multinom(k = k_multinom, options = options, A = A, b = b)
bayesfactor

