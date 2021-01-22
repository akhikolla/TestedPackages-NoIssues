# Test sparse version
if(T)
{
  # set.seed(123)
  N = 1000L
  d = 50L
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
  Xw = runif(N)
  clusterSizeUpperBound = sample(100L : 150L, length(centroidInd))


  sink("debug2.txt")
  rst = GMKMcharlie::sparseKMconstrained(
    X = sparsedat,
    # centroid = centroid,
    centroid = dat[centroidInd],
    Xw = Xw,
    clusterSizeUpperBound = clusterSizeUpperBound,
    minkP = 2,
    convergenceTail = 5L,
    tailConvergedRelaErr = 1e-5,
    maxIter = 100L,
    paraSortInplaceMerge = FALSE,
    maxCore = 7L,
    verbose = TRUE)
  rstTotalSS = sum(unlist(lapply(rst, function(x) sum(x$member2centroidDistance ^ 2))))
  rstClusterSizes = sort(unlist(lapply(rst, function(x) length(x$clusterMember))))
  sink()


  sink("debug2.txt")
  rst2 = GMKMcharlie::paraKMsizeConstrained(
    X = as.matrix(dat),
    centroid = as.matrix(dat[centroidInd]),
    Xw = Xw,
    clusterSizeUpperBound = clusterSizeUpperBound,
    minkP = 2,
    convergenceTail = 5L,
    tailConvergedRelaErr = 1e-5,
    maxIter = 100L,
    paraSortInplaceMerge = FALSE,
    maxCore = 7L,
    verbose = TRUE)
  sink()
  rstTotalSS2 = sum(unlist(lapply(rst, function(x) sum(x$member2centroidDistance ^ 2))))
  rstClusterSizes2 = sort(unlist(lapply(rst, function(x) length(x$clusterMember))))


  rstTotalSS / rstTotalSS2 - 1
  sum(abs(rstClusterSizes2 - rstClusterSizes))







  rst2 = kmeans(x = t(dat), centers = t(dat[, centroidInd]), iter.max = 100, algorithm = "Lloyd")
  rstTotalSS / rst2$tot.withinss - 1
  sum(abs(rstClusterSizes - sort(rst2$size)))
}




# Make final examples.
if(T)
{
  # ===========================================================================
  N = 5000L # Number of points.
  d = 500L # Dimensionality.
  K = 50L # Number of clusters.
  dat = matrix(rnorm(N * d) + runif(N * d), nrow = d)


  # Use kmeans++ initialization.
  centroidInd = GMKMcharlie::KMppIni(
    X, K, firstSelection = 1L, minkP = 2, stochastic = FALSE,
    seed = sample(1e9L, 1), maxCore = 2L, verbose = TRUE)


  centroid = dat[, centroidInd]


  # Each cluster should be no greater than N / K * 2.
  sizeConstraints = as.integer(rep(N / K * 2, K))
  system.time({rst = GMKMcharlie::KMconstrained(
    X = dat, centroid = centroid, clusterWeightLim = sizeConstraints,
    tailConvergedRelaErr = 1e-6, verbose = TRUE)})


  # Size upper bounds vary between N / K * 1.5 and N / K * 2
  sizeConstraints = as.integer(round(runif(K, N / K * 1.5, N / K * 2)))
  system.time({rst = GMKMcharlie::KMconstrained(
    X = dat, centroid = centroid, clusterWeightLim = sizeConstraints,
    tailConvergedRelaErr = 1e-6, verbose = TRUE)})


}

































