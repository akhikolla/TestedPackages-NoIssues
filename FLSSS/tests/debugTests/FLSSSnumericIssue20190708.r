

ds = sort(c(100, -50, -50, -100, -201.98, -532.31, -203.61, -2282.33, -2511.72, -73.75, -234.22, -494.73, -76.54, -297.9, -2392.27, -30.6, -98.33, -385.73, -132.18, -153.32, -114.36, -115, -256.88, -1157.27, 203.61, -217.8, -375, -289.06, -670.34, 55.2, -530.22, 378.21, -378.21, 670.34, 73.75, 256.88, 289.06, 217.8, 494.73, 2392.27, 76.54, 168.07, 129.83, 360.01, 115, 1157.27, -696.48, -30.19, 6030.57, 257.69, 97.37, -257.69, 530.22, 2282.33, 480.98, 51.33, 30.19, 234.22, 30.6, 153.32, 385.73, 132.18, 98.33, 114.36, 2511.72, -6030.57, -544.23, -78.66, 66.25, -66.25, 327.43, -327.43, -73.75, -232.91, 232.91, 73.75, 1097.5, 78.66, 544.23, -294.84, -203.81, -2282.33, 322.13, -322.13, 696.48, -1097.5, -55.2, -390.84, -234.92, -907.5))


# FLSSS:::z_findBound(len = 3, V = as.matrix(ds), target = 0, me = 0.001)
FLSSS:::z_findBound(len = 3, V = as.matrix(ds), target = 0, me = 0.001, initialLB = c(8,25,78)+1, initialUB = c(24,66,78)+1)


# lis = FLSSS::FLSSS(len = 3, v = ds * 1000, target = 0, ME = 0.001 * 1000, solutionNeed = 10)
# sink("debug.csv")
lis = FLSSS::FLSSS(len = 3, v = ds, target = 0, ME = 0.001, solutionNeed = 10, useBiSrchInFB = F)
# sink()
unlist(lapply(lis, function(x) sum(ds[x])))




library(FLSSS)

L = 5L # subset size

d = 5L # number of dimensions

N = 1000 # superset size

# generate a multidimensional vector/set
v = as.data.frame(matrix(abs(rnorm(d * N) * 1000), ncol = d))

# make an underlying solution
solution = sample(1L : N, L)

# the subset sum target of the underlying solution
target = as.numeric(colSums(v[solution, ]))

# bound the error as 5% of each dimension's target magnitude
ME = abs(target) * 0.05

# randomOrder = sample(1L:((nrow(v) - L) * L + 1L), (nrow(v) - L) * L + 1L)
# randomizeTargetOrder should be TRUE or FALSE or
# an index vector of size (N-L)*L+1

tmp = FLSSS::mFLSSSpar(maxCore = 7, L, v, target, ME, tlimit = 1000, solutionNeed = 50)




len = 3L
N = 10L
v = as.matrix(sort(rnorm(N)))
target = sum(sample(v, len))
me = 0.1 # error
rst = FLSSS:::z_findBound(len, v, target, me, initialLB = c(1, 2, 3), initialUB = c(8, 9, 10))
rst = rst[-1]
print(rst)


















