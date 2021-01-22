

N = 50
ndim = 5
len = 10


v = matrix(runif(N * ndim) * 1000, ncol = ndim)
v = apply(v, 2, function(x) sort(x))
sol = sort(sample(nrow(v), len))
target = colSums(v[sol, ])
me = abs(target) * 0.00001
tmp = FLSSS::mFLSSSpar(maxCore = 13, len = len, mV = v, mTarget = target, mME = me, solutionNeed = 1000L, tlimit = 60, avgThreadLoad = 5000)
for(i in 1:length(tmp)) tmp[[i]] = sort(tmp[[i]])


unlist(lapply(tmp, function(x) abs(colSums(v[x,]) - target) <= me))







