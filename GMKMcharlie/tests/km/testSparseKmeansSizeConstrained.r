

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
  Xw = numeric(0)
  rst = GMKMcharlie::paraSparseKMsizeConstrained(
    X = sparsedat,
    centroid = centroid,
    Xw = Xw,
    clusterSizeUpperBound = 50L,
    minkP = 2,
    convergenceTail = 5L,
    tailConvergedRelaErr = 1e-6,
    maxIter = 100L,
    paraSortInplaceMerge = FALSE,
    maxCore = 7L,
    verbose = TRUE)
  rstTotalSS = sum(unlist(lapply(rst, function(x) sum(x$member2centroidDistance ^ 2))))
  rstClusterSizes = sort(unlist(lapply(rst, function(x) length(x$clusterMember))))
  rst2 = kmeans(x = t(dat), centers = t(dat[, centroidInd]), iter.max = 100, algorithm = "Lloyd")
  rstTotalSS / rst2$tot.withinss - 1
  sum(abs(rstClusterSizes - sort(rst2$size)))
}





# Final examples to make.
if(T)
{
  # ===========================================================================
  # Play random numbers. See speed.
  # ===========================================================================
  N = 5000L # Number of points.
  d = 500L # Dimensionality.
  K = 50L # Number of clusters.


  # Create a data matrix, about 95% of which are zeros.
  dat = matrix(unlist(lapply(1L : N, function(x)
  {
    tmp = numeric(d)
    # Nonzero entries.
    Nnz = as.integer(max(1, d * runif(1, 0, 0.05)))
    tmp[sample(d, Nnz)] = runif(Nnz) + rnorm(Nnz)
    tmp
  })), nrow = d); gc()


  # Convert to sparse representation.
  # GMKMcharlie::d2s() is equivalent.
  sparsedat = apply(dat, 2, function(x)
  {
    nonz = which(x != 0)
    list(nonz, x[nonz])
  }); gc()


  centroidInd = sample(length(sparsedat), K)


  # Test speed using sparse representation.
  sparseCentroid = sparsedat[centroidInd]
  # Size upper bounds vary in [N / K * 1.5, N / K * 2]
  sizeConstraints = as.integer(round(runif(K, N / K * 1.5, N / K * 2)))
  system.time({sparseRst = GMKMcharlie::KMconstrainedSparse(
    X = sparsedat, d = d, centroid = sparseCentroid,
    clusterWeightLim = sizeConstraints,
    tailConvergedRelaErr = 1e-6,
    maxIter = 100, minkP = 2, maxCore = 2, verbose = TRUE)})




}


























