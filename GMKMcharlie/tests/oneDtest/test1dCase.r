

dat = seq(-12, 12, len = 1000)
w = 0.2 * dnorm(dat) + 0.3 * dnorm(dat, mean = -6, sd = 2) + 0.5 * dnorm(dat, mean = 5, sd = 3)
g = 3L
alpha = rep(1 / g, g)
mu = sample(as.numeric(dat), g)
s = rep(var(as.numeric(dat)) / g, g)
tmp = GMKMcharlie::paraGmm(X = t(dat), Xw = w, G = g, alpha = alpha, mu = t(mu), Sigma = t(s), maxIter = 10000, eigenRatioLim = 10)
plot(tmp$fitted)
tmp = GMKMcharlie::paraGmmCW(X = t(dat), Xw = w, G = g, alpha = alpha, mu = t(mu), Sigma = t(s), maxIter = 10000, eigenRatioLim = 10)
plot(tmp$fitted)
tmp = GMKMcharlie::paraGmmFJ(X = t(dat), Xw = w, G = g * 3, alpha = alpha, mu = t(mu), Sigma = t(s), maxIter = 10000, eigenRatioLim = 10)
plot(tmp$fitted)






















