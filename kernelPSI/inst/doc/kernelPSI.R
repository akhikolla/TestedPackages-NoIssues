## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE--------------------------------------------------------------
set.seed(14)

## -----------------------------------------------------------------------------
require("MASS")

n_kernels <- 10 # total number of kernels
size_kernels <- 5 # dimensionality of the data associated to each kernel

n <- 100 # sample size
rho <- 0.6 # correlation parameter

# intra-group correlation matrix
corr <- outer(seq_len(size_kernels), seq_len(size_kernels),
              function(i, j) return(rho^(abs(i-j))))

# design matrix
X <- replicate(n_kernels,
               mvrnorm(n, mu = rep(0, size_kernels), Sigma = corr),
               simplify = FALSE)

## -----------------------------------------------------------------------------
require("kernlab")
K <- replicate(n_kernels, rbfdot(sigma = 1 / size_kernels))  # full list of Gaussian kernels

# List of Gram matrices
Kmat <- sapply(seq_len(n_kernels),
               function(i) {kMatrix <- kernelMatrix(K[[i]], X[[i]]);
               return(as.kernelMatrix(kMatrix, center = TRUE))},
               simplify = FALSE)

## -----------------------------------------------------------------------------
m_kernels <- 3 # number of causal kernels
theta <- 0.1 # amplitude of size effect

Ksum <- Reduce(`+`, Kmat[seq_len(m_kernels)]) # sum kernel of the causal kernels
decompK <- eigen(Ksum) # eigenvalue decomposition of the sum kernel Ksum

Y <- as.matrix(theta * decompK$values[1] * decompK$vectors[, 1] + rnorm(n), ncol = 1) # response vector
Lmat <- kernelMatrix(new("vanillakernel"), Y) # linear response vector

## -----------------------------------------------------------------------------
require("kernelPSI")

candidate_kernels <- 3 # number of kernels for the fixed variant
selectFOHSIC <- FOHSIC(Kmat, Lmat, mKernels = candidate_kernels) # fixed variant
constraintFO <- forwardQ(Kmat, selectFOHSIC) # list of quadratic constraints modeling the selection event

## -----------------------------------------------------------------------------
selectAHSIC <- adaFOHSIC(Kmat, Lmat) # adaptive variant
adaFO <- adaQ(Kmat, selectAHSIC[["selection"]], selectAHSIC[["n"]]) # list of quadratic constraints for the adaptive selection method
adaS <- selectAHSIC$selection[seq_len(selectAHSIC$n)] # indices of selected kernels

## -----------------------------------------------------------------------------
n_replicates <- 5000 # number of replicates (statistical power and validity require a higher number of samples)
burn_in <- 1000 # number of burn-in iterations

# Fixed variant ------------------
# selected methods: 'ridge' for the kernel ridge regression prototype 
# and 'pca' for the kernel principal component regression prototype
kernelPSI(Y, K_select = Kmat[selectFOHSIC], constraintFO, method = c("ridge", "pca"), 
          n_replicates = n_replicates, burn_in = burn_in)

# Adaptive variant ------------------
# selected methods: 'hsic' for the unbiased HSIC estimator
kernelPSI(Y, K_select = Kmat[adaS], constraintFO, method = "hsic", 
          n_replicates = n_replicates, burn_in = burn_in)

