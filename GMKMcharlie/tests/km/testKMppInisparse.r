

# =============================================================================
N = 2000L
d = 3000L
X = matrix(rnorm(N * d) + 2, nrow = d)
# Fill many zeros in X:
X = apply(X, 2, function(x) {
  x[sort(sample(d, d * runif(1, 0.95, 0.99)))] = 0; x})
# Get the sparse version of X.
sparseX = GMKMcharlie::d2s(X)


K = 30L
seed = 123L
# Time cost of finding the centroids via dense representation.
# CRAN check allows only 2 threads. Increase `maxCore` for more speed.
system.time({kmppViaDense = GMKMcharlie::KMppIni(
  X, K, firstSelection = 1L, minkP = 2, stochastic = TRUE, seed = seed,
  maxCore = 2L)})


# Time cost of finding the initial centroids via sparse representation.
system.time({kmppViaSparse = GMKMcharlie::KMppIniSparse(
  sparseX, d, K, firstSelection = 1L, minkP = 2, stochastic = TRUE,
  seed = seed, maxCore = 2L)})


# Results should be identical.
sum(kmppViaSparse - kmppViaDense)












