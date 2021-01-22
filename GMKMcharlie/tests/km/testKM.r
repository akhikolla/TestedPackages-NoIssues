

if(T)
{
  set.seed(123)


  while(T) {
  # for(i in 1L : 100L)
  N = 5000L
  d = 30L
  K = 30L


  dat = matrix(unlist(lapply(1L : N, function(x)
  {
    tmp = numeric(d)
    Nnz = as.integer(max(1, d * runif(1, 0, 0.25)))
    # Nnz = 2L
    tmp[sample(d, Nnz)] = runif(Nnz) + rnorm(Nnz)
    tmp
  })), nrow = d); gc()


  if(F)
  {
    dat = apply(dat, 2, function(x) x / sum(x ^ 2) ^ 0.5)
  }; gc()


  sparsedat = apply(dat, 2, function(x)
  {
    nonz = which(x != 0)
    list(nonz, x[nonz])
  }); gc()


  centroidInd = sample(length(sparsedat), K)
  centroid = sparsedat[centroidInd]


  # save.image()


  # sink("debug.txt")
  system.time({sparseRst = GMKMcharlie::KMsparse(
    X = sparsedat,
    d = d,
    centroid = centroid,
    minkP = "max",
    maxIter = 100L,
    maxCore = 1L,
    verbose = TRUE)})
  # sink()


  sparseRstTotalSS = sum(unlist(lapply(sparseRst, function(x) sum(x$member2centroidDistance ^ 2))))
  sparseRstClusterSizes = sort(unlist(lapply(sparseRst, function(x) length(x$clusterMember))))


  # sink("debug2.txt")
  system.time({rst = GMKMcharlie::KM(
    X = as.matrix(dat),
    centroid = as.matrix(dat[, centroidInd]),
    minkP = "max",
    maxIter = 100L,
    maxCore = 1L,
    verbose = TRUE)})
  # sink()


  rstTotalSS = sum(unlist(lapply(rst, function(x) sum(x$member2centroidDistance ^ 2))))
  rstClusterSizes = sort(unlist(lapply(rst, function(x) length(x$clusterMember))))

  tmp = kmeans(t(dat), t(dat[,centroidInd]), algorithm = "Lloyd")


  if(abs(sparseRstTotalSS / rstTotalSS - 1) > 1e-5) break


  }







  system.time({rstR = kmeans(x = t(dat), centers = t(dat[, centroidInd]), iter.max = 100, algorithm = "Lloyd", trace = TRUE)})
  system.time({rstR = testStatsKM(x = t(dat), centroid = t(dat[, centroidInd]), maxIter = 100L)})
  rstTotalSS / rstR$tot.withinss - 1
  sum(abs(rstClusterSizes - sort(rstR$size)))
}




# Test kmeans with cluster weight constraints.
if(T)
{
  set.seed(123)


  # while(T) {
  N = 30000L
  d = 5000L
  K = 100L
  dat = matrix(rnorm(N * d) + runif(N * d), ncol = d)
  dat = as.data.frame(apply(dat, 1, function(x)
  {
    r = runif(1, 0.9, 1)
    l = length(x)
    x[sample(l, min(d - 1, l * r))] = 0
    x
  }))
  sparsedat = lapply(dat, function(x)
  {
    nonz = which(x != 0)
    list(nonz, x[nonz])
  })
  centroidInd = sample(length(sparsedat), K)
  centroid = sparsedat[centroidInd]


  # save.image()


  # sink("debug.txt")
  dataPointWeights = abs(rnorm(length(sparsedat)))
  dataPointWeights = dataPointWeights / sum(dataPointWeights) * length(sparsedat)
  clusterWeightLims = sample(1000 : 2000, length(centroid), replace = TRUE)
  system.time({sparseRst = GMKMcharlie::sparseKMconstrained(
    X = sparsedat,
    d = d,
    centroid = centroid,
    # Xw = dataPointWeights,
    clusterWeightLim = rep(350, length(centroid)),
    minkP = 2,
    convergenceTail = 5L,
    tailConvergedRelaErr = 1e-10,
    maxIter = 1000L,
    maxCore = 7L,
    verbose = TRUE)})
  # sink()


  sparseRstTotalSS = sum(unlist(lapply(sparseRst, function(x) sum(x$member2centroidDistance ^ 2))))
  sparseRstClusterSizes = sort(unlist(lapply(sparseRst, function(x) length(x$clusterMember))))


  # sink("debug2.txt")
  system.time({rst = GMKMcharlie::KMconstrained(
    X = as.matrix(dat),
    centroid = as.matrix(dat[centroidInd]),
    clusterWeightLim = rep(350, length(centroidInd)),
    minkP = 2,
    convergenceTail = 5L,
    tailConvergedRelaErr = 1e-10,
    maxIter = 1000L,
    maxCore = 7L,
    verbose = TRUE)})
  # sink()


  rstTotalSS = sum(unlist(lapply(rst, function(x) sum(x$member2centroidDistance ^ 2))))
  rstClusterSizes = sort(unlist(lapply(rst, function(x) length(x$clusterMember))))
  # if(abs(sparseRstTotalSS / rstTotalSS - 1) > 1e-5) break
  # }







  system.time({rstR = kmeans(x = t(dat), centers = t(dat[, centroidInd]), iter.max = 100, algorithm = "Lloyd", trace = TRUE)})
  system.time({rstR = testStatsKM(x = t(dat), centroid = t(dat[, centroidInd]), maxIter = 100L)})
  rstTotalSS / rstR$tot.withinss - 1
  sum(abs(rstClusterSizes - sort(rstR$size)))
}




# Test kmeans with weights.
if(T)
{
  # seed = sample(1e9L, 1)
  # seed
  # set.seed(seed)
  set.seed(334823933)
  K = 10L
  X = as.matrix(expand.grid(seq(-10, 10, len = 100), seq(-10, 10, len = 100)))
  dimnames(X) = NULL
  para = lapply(1 : K, function(x)
  {
    mu = runif(2, -9, 9)
    s = runif(2, 1, 3)
    rho = runif(1, -0.8, 0.8)
    s = s %*% t(s)
    s[2:3] = s[2:3] * rho
    list(m = mu, s = s)
  })
  w = rowSums(as.data.frame(lapply(para, function(x)
  {
    mvtnorm::dmvnorm(X, mean = x$m, sigma = x$s)
  })))
  image(matrix(w, ncol = 100))
  names(w) = NULL


  sparsedat = apply(X, 1, function(x) data.frame(c(1L, 2L), x))
  centroidInd = sample(length(sparsedat), 10)
  rst = GMKMcharlie::paraSparseKMclassic(X = sparsedat, centroid = sparsedat[centroidInd], Xw = w, maxIter = 100, maxCore = 7L, verbose = TRUE)
}




if(T)
{
  set.seed(123)
  N = 1000L
  d = 5L
  K = 20L
  dat = matrix(rnorm(N * d) + runif(N * d), ncol = d)
  dat = as.data.frame(apply(dat, 1, function(x)
  {
    r = runif(1, 0.5, 1)
    l = length(x)
    x[sample(l, min(d - 1, l * r))] = 0
    x
  }))
  sparsedat = lapply(dat, function(x)
  {
    nonz = which(x != 0)
    list(which(x != 0), x[nonz])
  })
  centroidInd = sample(length(sparsedat), K)
  centroid = sparsedat[centroidInd]
  # Xw = sample(5, length(sparsedat), replace = T)
  Xw = numeric(0)
  rst = GMKMcharlie::paraSparseKMclassic(X = sparsedat, centroid = centroid, Xw = Xw, maxIter = 100, maxCore = 7L, verbose = TRUE)
  rstTotalSS = sum(unlist(lapply(rst, function(x) sum(x$member2centroidDistance ^ 2))))
  rstClusterSizes = sort(unlist(lapply(rst, function(x) length(x$clusterMember))))


  # tmp = mapply(function(x, y) rep(x, y), dat, Xw, SIMPLIFY = F)
  # tmp = t(matrix(unlist(tmp), nrow = d))
  tmp = t(as.data.frame(tmp))
  rst2 = kmeans(x = tmp, centers = t(dat[centroidInd]), iter.max = 100, algorithm = "Lloyd")
  rstTotalSS / rst2$tot.withinss - 1
  sum(abs(rstClusterSizes - sort(rst2$size)))
}




# Test dense Kmeans
if(T)
{
  N = 1000L
  d = 10L
  K = 20L
  dat = t(matrix(rnorm(N * d) + runif(N * d), ncol = d))
  centroidInd = sample(ncol(dat), K)
  centroid = dat[, centroidInd]
  rst = GMKMcharlie::paraKMclassic(X = dat, centroid = centroid, Xw = numeric(0), maxIter = 100, maxCore = 7L, verbose = TRUE)
  rst2 = kmeans(x = t(dat), centers = t(centroid), iter.max = 100, algorithm = "Lloyd")
  sum(unlist(lapply(rst, function(x) sum(x$member2centroidDistance ^ 2)))) / rst2$tot.withinss - 1


  K = 20L
  dat = t(iris[1:4])
  centroidInd = sample(ncol(dat), K)
  centroid = dat[, centroidInd]
  rst = GMKMcharlie::paraKMclassic(X = dat, centroid = centroid, Xw = numeric(0), maxIter = 100, maxCore = 7L, verbose = TRUE)
  rst2 = kmeans(x = t(dat), centers = t(centroid), iter.max = 100, algorithm = "Lloyd")
  sum(unlist(lapply(rst, function(x) sum(x$member2centroidDistance ^ 2)))) / rst2$tot.withinss - 1
}




# Final examples to make.
if(T)
{
  # ===========================================================================
  # Play random numbers. See speed.
  # ===========================================================================
  N = 10000L # Number of points.
  d = 500L # Dimensionality.
  K = 100L # Number of clusters.
  dat = matrix(rnorm(N * d) + runif(N * d), nrow = d)


  # Use kmeans++ initialization.
  centroidInd = GMKMcharlie::KMppIni(
    X, K, firstSelection = 1L, minkP = 2, stochastic = FALSE,
    seed = sample(1e9L, 1), maxCore = 2L, verbose = TRUE)


  centroid = dat[, centroidInd]


  # Euclidean.
  system.time({rst = GMKMcharlie::KM(
    X = dat, centroid = centroid, maxIter = 100,
    minkP = 2, maxCore = 2, verbose = TRUE)})


  # Cosine dissimilarity.
  dat = apply(dat, 2, function(x) x / sum(x ^ 2) ^ 0.5)
  centroid = dat[, centroidInd]
  system.time({rst2 = GMKMcharlie::KM(
    X = dat, centroid = centroid, maxIter = 100,
    minkP = "cosine", maxCore = 2, verbose = TRUE)})




  # ===========================================================================
  # Test against R's inbuilt km()
  # ===========================================================================
  dat = t(iris[1:4])
  dimnames(dat) = NULL


  # Use kmeans++ initialization.
  centroidInd = GMKMcharlie::KMppIni(
    X = dat, K = 3, firstSelection = 1L, minkP = 2, stochastic = FALSE,
    seed = sample(1e9L, 1), maxCore = 2L, verbose = TRUE)
  centroid = dat[, centroidInd]


  rst = GMKMcharlie::KM(X = dat, centroid = centroid, maxIter = 100,
                        minkP = 2, maxCore = 2, verbose = TRUE)
  rst = lapply(rst, function(x) sort(x$clusterMember))


  rst2 = kmeans(x = t(dat), centers = t(centroid), algorithm = "Lloyd")
  rst2 = aggregate(list(1L : length(rst2$cluster)),
                   list(rst2$cluster), function(x) sort(x))[[2]]


  setdiff(rst, rst2)


}































