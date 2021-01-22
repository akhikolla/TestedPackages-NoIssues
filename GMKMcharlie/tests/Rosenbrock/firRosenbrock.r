

# =============================================================================
# set.seed(12345)
f <- function(x, y) { (1 - x) ^ 2 + 100 * (y - x ^ 2) ^ 2 }
tmp = expand.grid(seq(-1.5, 2, len = 50), seq(-0.5, 3, len = 50))
x = tmp[[1]]; y = tmp[[2]]
# x = runif(1000, -1.5, 2); y = runif(1000, -0.5, 3)
z = f(x, y)
x = x / sd(x)
y = y / sd(y)
z = z / sd(z)
dat = rbind(x, y)
dimnames(dat) = NULL


# sink("debug.txt")
fittedGM = GMKMcharlie::paraGmm(
  dat, G = 5, Xw = z, maxCore = 7, maxIter = 1000L,
  eigenRatioLim = 0, convergenceEPS = 1e-5, alphaEPS = 1e-16,
  tlimit = 3600, verbose = 10)
# sink()
rgl::plot3d(x = x, y = y, z = fittedGM$fitted)




fittedGM = GMKMcharlie::paraGmmCW(
  dat, G = 6, Xw = z, maxCore = 7, maxIter = 1000,
  eigenRatioLim = 0, convergenceEPS = 1e-5, alphaEPS = 1e-16,
  tlimit = 3600, verbose = 10)
# sink()
rgl::plot3d(x = x, y = y, z = fittedGM$fitted)




fittedGM = GMKMcharlie::paraGmmFJ(
  dat, G = 100, Gmin = 2, Xw = z, maxCore = 7, maxIter = 1000,
  eigenRatioLim = 0, convergenceEPS = 1e-5, alphaEPS = 1e-16,
  tlimit = 3600, verbose = 10)
# sink()
rgl::plot3d(x = x, y = y, z = fittedGM$fitted, type = "h", col = "red")





# Iris data
fittedGM = GMKMcharlie::paraGmmCW(
  t(iris[1:4]), G = 3, Xw = numeric(0), maxCore = 7, maxIter = 1000,
  eigenRatioLim = 0, convergenceEPS = 1e-5, alphaEPS = 1e-16,
  tlimit = 3600, verbose = 10)


fittedGM = GMKMcharlie::paraGmm(
  t(iris[1:4]), G = 3, Xw = numeric(0), maxCore = 7, maxIter = 1000,
  eigenRatioLim = 0, convergenceEPS = 1e-5, alphaEPS = 1e-16,
  tlimit = 3600, verbose = 10)


fittedGM = GMKMcharlie::paraGmmFJ(
  t(iris[1:4]), G = 5, Gmin = 2, Xw = numeric(0), maxCore = 1, maxIter = 1000,
  eigenRatioLim = 0, convergenceEPS = 1e-5, alphaEPS = 1e-16,
  tlimit = 3600, verbose = 10)










