

# =============================================================================
N = 30000L
d = 300L
X = matrix(rnorm(N * d) + 2, nrow = d)
K = 30L
kmppSt = GMKMcharlie::KMppIni(X, K, firstSelection = 1L, minkP = 2,
                 stochastic = TRUE, seed = sample(1e9L, 1), maxCore = 2L)
kmppDt = GMKMcharlie::KMppIni(X, K, firstSelection = 1L, minkP = 2,
                 stochastic = FALSE, maxCore = 2L)
str(kmppSt)
str(kmppDt)














# N = 30000L
# d = 300L
# X = matrix(rnorm(N * d) + 2, nrow = d)
# zeroR = runif(1, 0.8, 0.99)
# X[sample(N, N * zeroR)] = 0
# X = apply(X, 2, function(x) {if(sum(x) == 0) x[sample(d, 1)] = rnorm(1); x})
# system.time({Xsparse = d2s(X)})
# # system.time({Xsparse2 = apply(X, 2, function(x)
# # {
# #   nz = which(abs(x - 0) >= 1e-16)
# #   data.frame(nz, x[nz])
# # })})
# # X = s2d(Xsparse, d)
# K = 35L
# seed = sample(1e9L, 1)
# kmpp = KMppIni(X, K, firstSelection = 1L, minkP = 2,
#                stochastic = TRUE, seed = seed, maxCore = 2L)
# kmppSparse = KMppIniSparse(Xsparse, d, K, firstSelection = 1L, minkP = 2,
#                      stochastic = TRUE, seed = seed, maxCore = 7L)
# sum(abs(kmpp - kmppSparse))















