

# Function hub.
if(T)
{
  # x is column-major data matrix
  computeDensityMat <- function(x, glist, withAlpha = F)
  {
    a = glist$a
    m = glist$mu
    s = glist$sigma
    if(is.matrix(m)) d = nrow(m)
    else d = length(m)


    rst = as.matrix(as.data.frame(mapply(function(weight, mean, covar)
    {
      tmp =  mvtnorm::dmvnorm(x, mean = mean, sigma = matrix(covar, ncol = d))
      if(withAlpha) tmp * weight
      else tmp
    }, a, as.data.frame(m), as.data.frame(s), SIMPLIFY = F)))
    names(rst) = NULL
    rst
  }


  weightDenMat <- function(denMat, alpha)
  {
    rst = as.matrix(as.data.frame(mapply(function(x, y) x * y, as.data.frame(denMat), alpha, SIMPLIFY = F)))
    dimnames(rst) = NULL
    rst
  }


  if(F)
  {
    FJgmmOPT01 <- function(X, glist, kmin = 2L, eps = 1e-4) # Examined
    {
      thetaT_1 = glist
      thetaT_1loss = 1e308
      minLoss = 1e308
      Xsize = nrow(X)
      denMat = computeDensityMat(X, glist, withAlpha = T)
      t = 1L
      d = nrow(glist$mu)
      N = d + d * (d + 1) / 2
      while(T)
      {
        while(T)
        {
          cat("t =", t, "")
          t = t + 1L # Useless
          m = 1L
          while(m <= length(glist$a))
          {
            Wmat = denMat / rowSums(denMat)
            aArray = apply(Wmat, 2, function(x)
            {
              max(sum(x) - N / 2, 0)
            })
            glist$a[m] = aArray[m] / sum(aArray)
            glist$a = glist$a / sum(glist$a)


            if(glist$a[m] > 0)
            {
              pointw = Wmat[, m] / sum(Wmat[, m])
              glist$mu[, m] = colSums(X * pointw)
              centered = t(X) - glist$mu[, m]
              glist$sigma[, m] = as.numeric(centered %*% (t(centered) * pointw))
              m = m + 1L
            }
            else
            {
              glist$a = glist$a[-m]
              glist$mu = glist$mu[, -m]
              glist$sigma = glist$sigma[, -m]
            }


            denMat = computeDensityMat(X, glist, withAlpha = T)
          }


          thetaT = glist # Useless, just notation
          knz = length(glist$a)
          thetaTloss = N / 2 * sum(log(thetaT$a * Xsize / 12)) + knz / 2 * log(Xsize / 12) + (knz * N + knz) / 2 - sum(log(rowSums(denMat)))


          cat("knz =", knz, "")
          cat("prior loss =", thetaT_1loss, "")
          cat("loss =", thetaTloss, "")
          cat("loss relative difference =", thetaTloss / thetaT_1loss - 1, "\n")


          if(thetaT_1loss - thetaTloss < eps * abs(thetaT_1loss)) break
          thetaT_1loss = thetaTloss
        }
        print("Component-wise GMM converged")


        if(thetaTloss <= minLoss)
        {
          minLoss = thetaTloss
          bestTheta = glist
        }


        whichMin = which.min(glist$a)
        glist$a = glist$a[-whichMin]
        glist$mu = glist$mu[, -whichMin]
        glist$sigma = glist$sigma[, -whichMin]


        if(length(glist$a) >= kmin) denMat = computeDensityMat(X, glist, withAlpha = T)
        else break
      }
      bestTheta
    }
  }


  # Loss comparison has no absolute value constraint..
  FJgmmNoAbs <- function(X, glist, kmin = 2L, eps = 1e-5)
  {
    thetaT_1 = glist
    thetaT_1loss = 1e308
    minLoss = 1e308
    Xsize = nrow(X)
    denMat = computeDensityMat(X, glist, withAlpha = F)
    t = 1L
    d = nrow(glist$mu)
    N = d + d * (d + 1) / 2
    while(T)
    {
      while(T)
      {
        # cat("t =", t, "\n")
        t = t + 1L # Useless
        m = 1L
        while(m <= length(glist$a))
        {
          denMatWithAlpha = weightDenMat(denMat, glist$a)
          denMatWithAlphaRowSums = rowSums(denMatWithAlpha)


          # cat("rowSums = \n")
          # print(denMatWithAlphaRowSums)


          Wmat = denMatWithAlpha / denMatWithAlphaRowSums
          # cat("\n==============================\n")
          # print(Wmat)
          # cat("==============================\n")


          aArray = apply(Wmat, 2, function(x)
          {
            max(sum(x) - N / 2, 0)
          })
          # sumArray = sum(aArray)
          # if(sumArray == 0)
          glist$a[m] = aArray[m] / sum(aArray)
          glist$a = glist$a / sum(glist$a)


          cat("All alphas:", glist$a, "\n\n")
          if(glist$a[m] > 0)
          {
            pointw = Wmat[, m] / sum(Wmat[, m])
            glist$mu[, m] = colSums(X * pointw)
            centered = t(X) - glist$mu[, m]
            glist$sigma[, m] = as.numeric(centered %*% (t(centered) * pointw))
            m = m + 1L
          }
          else
          {
            # knz = knz - 1L
            glist$a = glist$a[-m]
            glist$mu = glist$mu[, -m]
            glist$sigma = glist$sigma[, -m]
          }


          # cat(glist$a, "\n")


          denMat = computeDensityMat(X, glist, withAlpha = F)
        }


        thetaT = glist # Useless, just notation
        denMatWithAlphaRowSums = rowSums(weightDenMat(denMat, glist$a))
        knz = length(glist$a)
        thetaTloss = N / 2 * sum(log(thetaT$a[thetaT$a > 0] * Xsize / 12)) + knz / 2 * log(Xsize / 12) + (knz * N + knz) / 2 - sum(log(denMatWithAlphaRowSums))


        cat("prior loss =", thetaT_1loss, "")
        cat("loss =", thetaTloss, "\n")


        # if(thetaT_1loss >= thetaTloss & abs(thetaT_1loss - thetaTloss) < eps * abs(thetaT_1loss)) break
        if(thetaT_1loss - thetaTloss < eps * abs(thetaT_1loss)) break
        # if(abs(thetaT_1loss - thetaTloss) < eps * abs(thetaT_1loss)) break
        thetaT_1loss = thetaTloss
      }
      print("Component-wise GMM converged")


      if(thetaTloss <= minLoss)
      {
        # cat("thetaTloss <= minLoss\n")
        minLoss = thetaTloss
        bestTheta = glist
      }


      whichMin = which.min(glist$a)
      glist$a = glist$a[-whichMin]
      glist$a = glist$a / sum(glist$a)
      glist$mu = glist$mu[, -whichMin]
      glist$sigma = glist$sigma[, -whichMin]


      if(length(glist$a) >= kmin)
      {
        # denMat = computeDensityMat(X, glist, withAlpha = F)
        denMat = denMat[, -whichMin]
      }
      else break
    }




    denWithAlpha = computeDensityMat(X, bestTheta, withAlpha = T)
    tmp = apply(denWithAlpha, 1, function(x) which.max(x))
    clusterMember = aggregate(list(1L : nrow(X)), list(tmp), function(x) x)[[2]]
    for(i in 1L : length(clusterMember)) clusterMember[[i]] = sort(clusterMember[[i]])


    list(bestTheta = bestTheta, fitted = apply(computeDensityMat(X, bestTheta, withAlpha = T), 1, function(x) sum(x)), minLoss = minLoss, clusterMember = clusterMember)
  }


  FJgmm <- function(X, glist, kmin = 2L, eps = 1e-5)
  {
    thetaT_1 = glist
    thetaT_1loss = 1e308
    minLoss = 1e308
    Xsize = nrow(X)
    denMat = computeDensityMat(X, glist, withAlpha = F)
    t = 1L
    d = nrow(glist$mu)
    N = d + d * (d + 1) / 2
    while(T)
    {
      while(T)
      {
        # cat("t =", t, "\n")
        t = t + 1L # Useless
        m = 1L
        while(m <= length(glist$a))
        {
          denMatWithAlpha = weightDenMat(denMat, glist$a)
          denMatWithAlphaRowSums = rowSums(denMatWithAlpha)


          # cat("rowSums = \n")
          # print(denMatWithAlphaRowSums)


          Wmat = denMatWithAlpha / denMatWithAlphaRowSums
          # cat("\n==============================\n")
          # print(Wmat)
          # cat("==============================\n")


          aArray = apply(Wmat, 2, function(x)
          {
            max(sum(x) - N / 2, 0)
          })
          # sumArray = sum(aArray)
          # if(sumArray == 0)
          glist$a[m] = aArray[m] / sum(aArray)
          glist$a = glist$a / sum(glist$a)


          # cat("All alphas:", glist$a, "\n\n")
          if(glist$a[m] > 0)
          {
            pointw = Wmat[, m] / sum(Wmat[, m])
            glist$mu[, m] = colSums(X * pointw)
            centered = t(X) - glist$mu[, m]
            glist$sigma[, m] = as.numeric(centered %*% (t(centered) * pointw))
            m = m + 1L
          }
          else
          {
            # knz = knz - 1L
            glist$a = glist$a[-m]
            glist$mu = glist$mu[, -m]
            glist$sigma = glist$sigma[, -m]
          }


          # cat(glist$a, "\n")


          denMat = computeDensityMat(X, glist, withAlpha = F)
        }


        thetaT = glist # Useless, just notation
        denMatWithAlphaRowSums = rowSums(weightDenMat(denMat, glist$a))
        knz = length(glist$a)
        thetaTloss = N / 2 * sum(log(thetaT$a[thetaT$a > 0] * Xsize / 12)) + knz / 2 * log(Xsize / 12) + (knz * N + knz) / 2 - sum(log(denMatWithAlphaRowSums))


        cat("prior loss =", thetaT_1loss, "")
        cat("loss =", thetaTloss, "\n")


        # if(thetaT_1loss >= thetaTloss & abs(thetaT_1loss - thetaTloss) < eps * abs(thetaT_1loss)) break
        # if(thetaT_1loss - thetaTloss < eps * abs(thetaT_1loss)) break
        if(abs(thetaT_1loss - thetaTloss) < eps * abs(thetaT_1loss)) break
        thetaT_1loss = thetaTloss
      }
      print("Component-wise GMM converged")


      if(thetaTloss <= minLoss)
      {
        # cat("thetaTloss <= minLoss\n")
        minLoss = thetaTloss
        bestTheta = glist
      }


      whichMin = which.min(glist$a)
      glist$a = glist$a[-whichMin]
      glist$a = glist$a / sum(glist$a)
      glist$mu = glist$mu[, -whichMin]
      glist$sigma = glist$sigma[, -whichMin]


      if(length(glist$a) >= kmin)
      {
        # denMat = computeDensityMat(X, glist, withAlpha = F)
        denMat = denMat[, -whichMin]
      }
      else break
    }




    denWithAlpha = computeDensityMat(X, bestTheta, withAlpha = T)
    tmp = apply(denWithAlpha, 1, function(x) which.max(x))
    clusterMember = aggregate(list(1L : nrow(X)), list(tmp), function(x) x)[[2]]
    for(i in 1L : length(clusterMember)) clusterMember[[i]] = sort(clusterMember[[i]])


    list(bestTheta = bestTheta, fitted = apply(computeDensityMat(X, bestTheta, withAlpha = T), 1, function(x) sum(x)), minLoss = minLoss, clusterMember = clusterMember)
  }


  FJgmmOriginal <- function(X, glist, kmin = 2L, eps = 1e-5)
  {
    thetaT_1 = glist
    thetaT_1loss = 1e308
    minLoss = 1e308
    Xsize = nrow(X)
    denMat = computeDensityMat(X, glist, withAlpha = F)
    t = 1L
    d = nrow(glist$mu)
    N = d + d * (d + 1) / 2
    kmax = sum(glist$a > 0)
    knz = kmax


    while(knz >= kmin)
    {
      while(knz >= kmin)
      # while(T)
      {
        # cat("t =", t, "\n")
        t = t + 1L # Useless
        # m = 1L


        for(m in 1L : kmax)
        {
          if(knz < kmin) break


          denMatWithAlpha = weightDenMat(denMat, glist$a)
          # denMatWithAlpha = computeDensityMat(X, glist, withAlpha = T)
          denMatWithAlphaRowSums = rowSums(denMatWithAlpha)


          Wmat = denMatWithAlpha / denMatWithAlphaRowSums


          aArray = apply(Wmat, 2, function(x)
          {
            max(sum(x) - N / 2, 0)
          })


          sumArray = sum(aArray)
          if(sumArray == 0) glist$a = numeric(length(glist$a))
          else
          {
            glist$a[m] = aArray[m] / sumArray
            glist$a = glist$a / sum(glist$a)
          }


          # cat("glist$a = ", glist$a, "\n")


          if(glist$a[m] > 0)
          {
            pointw = Wmat[, m] / sum(Wmat[, m])
            glist$mu[, m] = colSums(X * pointw)
            centered = t(X) - glist$mu[, m]
            glist$sigma[, m] = as.numeric(centered %*% (t(centered) * pointw))
            denMat = computeDensityMat(X, glist, withAlpha = F)
          }
          else
          {
            # knz = knz - 1L
            # knz = sum(glist$a > 0)
          }
          knz = sum(glist$a > 0)


          # cat(glist$a, "\n")


        }


        thetaT = glist # Useless, just notation
        denMatWithAlphaRowSums = rowSums(computeDensityMat(X, thetaT, withAlpha = T))
        # denMatWithAlphaRowSums = rowSums(weightDenMat(denMat, glist$a))
        # knz = length(glist$a)
        knz = sum(glist$a > 0)
        thetaTloss = N / 2 * sum(log(thetaT$a[thetaT$a > 0] * Xsize / 12)) + knz / 2 * log(Xsize / 12) + (knz * N + knz) / 2 - sum(log(denMatWithAlphaRowSums))
        # cat("knz =", knz, "\n")
        # cat("loss =", thetaTloss, "\n")
        cat("prior loss =", thetaT_1loss, "")
        cat("loss =", thetaTloss, "\n")


        # if(thetaT_1loss >= thetaTloss & abs(thetaT_1loss - thetaTloss) < eps * abs(thetaT_1loss)) break
        # if(thetaT_1loss - thetaTloss < eps * abs(thetaT_1loss)) break
        if(abs(thetaT_1loss - thetaTloss) < eps * abs(thetaT_1loss)) break
        thetaT_1loss = thetaTloss
      }


      print("Component-wise GMM converged")


      if(thetaTloss <= minLoss)
      {
        # cat("thetaTloss <= minLoss\n")
        minLoss = thetaTloss
        bestTheta = glist
      }


      tmp = glist$a
      tmp[tmp == 0] = 1e300
      whichMin = which.min(tmp)
      glist$a[whichMin] = 0
      # cat("knz =", knz, " --- ")
      # knz = knz - 1L
      knz = sum(glist$a > 0)
      # cat("knz =", knz, " --- ")
      glist$a = glist$a / sum(glist$a)
      # denMat = computeDensityMat(X, glist, withAlpha = F)
    }


    tmp = which(bestTheta$a > 0)
    bestTheta$a = bestTheta$a[tmp]
    bestTheta$mu = bestTheta$mu[, tmp]
    bestTheta$sigma = bestTheta$sigma[, tmp]


    denWithAlpha = computeDensityMat(X, bestTheta, withAlpha = T)
    tmp = apply(denWithAlpha, 1, function(x) which.max(x))
    clusterMember = aggregate(list(1L : nrow(X)), list(tmp), function(x) x)[[2]]
    for(i in 1L : length(clusterMember)) clusterMember[[i]] = sort(clusterMember[[i]])


    list(bestTheta = bestTheta, fitted = apply(computeDensityMat(X, bestTheta, withAlpha = T), 1, function(x) sum(x)), minLoss = minLoss, clusterMember = clusterMember)
  }


  # From FJ's matlab code. The only difference is a slight difference in the order of update in a loop. There is no errors in principle against FJ's algorithm.
  FJgmmMatlab <- function(X, glist, kmin = 2L, eps = 1e-4)
  {
    thetaT_1 = glist
    thetaT_1loss = 1e308
    loglike_1 = 1e308
    minLoss = 1e308
    Xsize = nrow(X)
    denMat = computeDensityMat(X, glist, withAlpha = F)
    t = 1L
    d = nrow(glist$mu)
    N = d + d * (d + 1) / 2
    kmax = sum(glist$a > 0)
    knz = kmax


    thetaTlossSeq = list()
    globalCounter = 0L
    while(knz >= kmin)
    {
      while(knz >= kmin)
        # while(T)
      {
        # cat("t =", t, "\n")
        t = t + 1L # Useless
        # m = 1L


        for(m in 1L : kmax)
        {
          globalCounter = globalCounter + 1L


          if(knz < kmin) break
          denMatWithAlpha = weightDenMat(denMat, glist$a)
          # denMatWithAlpha = computeDensityMat(X, glist, withAlpha = T)
          denMatWithAlphaRowSums = rowSums(denMatWithAlpha)


          Wmat = denMatWithAlpha / denMatWithAlphaRowSums


          if(F)
          {
            aArray = apply(Wmat, 2, function(x)
            {
              max(sum(x) - N / 2, 0)
            })


            sumArray = sum(aArray)
            if(sumArray == 0) glist$a = numeric(length(glist$a))
            else
            {
              glist$a[m] = aArray[m] / sumArray
              glist$a = glist$a / sum(glist$a)
            }
          }


          # ..
          glist$a[m] = max(sum(Wmat[, m]) - N / 2, 0) / Xsize
          glist$a = glist$a / sum(glist$a)

          # cat("glist$a = ", glist$a, "\n")


          if(glist$a[m] > 0)
          {
            pointw = Wmat[, m] / sum(Wmat[, m])
            glist$mu[, m] = colSums(X * pointw)
            centered = t(X) - glist$mu[, m]
            glist$sigma[, m] = as.numeric(centered %*% (t(centered) * pointw))
            # if(globalCounter <= 250 & globalCounter >= 185)
            # {
            #   print(glist$sigma)
            #   print("======================================")
            # }
            denMat = computeDensityMat(X, glist, withAlpha = F)
          }
          else
          {
            # knz = knz - 1L
            # knz = sum(glist$a > 0)
            # cat("knz =", sum(glist$a > 0), "  ")
          }
          knz = sum(glist$a > 0)


          # cat(glist$a, "\n")


        }


        thetaT = glist # Useless, just notation
        denMatWithAlphaRowSums = rowSums(computeDensityMat(X, thetaT, withAlpha = T))
        # cat(denMatWithAlphaRowSums, "\n")
        # denMatWithAlphaRowSums = rowSums(weightDenMat(denMat, glist$a))
        # knz = length(glist$a)
        loglike = sum(log(denMatWithAlphaRowSums))
        knz = sum(glist$a > 0)
        cat("knz =", knz, " \n")
        # thetaTloss = N / 2 * sum(log(thetaT$a[thetaT$a > 0] * Xsize / 12)) + knz / 2 * log(Xsize / 12) + (knz * N + knz) / 2 - loglike
        thetaTloss = -loglike + N / 2 * sum(log(thetaT$a[thetaT$a > 0])) + (N / 2 + 0.5) * knz * log(Xsize)
        thetaTlossSeq[[length(thetaTlossSeq) + 1]] = thetaTloss


        # cat("knz =", knz, "\n")
        # cat("loss =", thetaTloss, "\n")
        cat("Prior neg loglike =", -loglike_1, "")
        cat("Neg loglike =", -loglike, "\n")
        # cat("loglikePrint =", loglike, "\n")
        # cat("dlength =", thetaTloss, "\n")


        # if(thetaT_1loss >= thetaTloss & abs(thetaT_1loss - thetaTloss) < eps * abs(thetaT_1loss)) break
        # if(thetaT_1loss - thetaTloss < eps * abs(thetaT_1loss)) break
        # if(abs(thetaT_1loss - thetaTloss) < eps * abs(thetaT_1loss)) break
        # if(abs(loglike_1 - ))
        if(abs(loglike / loglike_1 - 1) < eps) break
        loglike_1 = loglike
      }


      print("Component-wise GMM converged")


      # if(thetaTloss <= minLoss)
      if(thetaTlossSeq[[length(thetaTlossSeq)]] < minLoss)
      {
        # cat("thetaTloss <= minLoss\n")
        # minLoss = thetaTloss
        minLoss = thetaTlossSeq[[length(thetaTlossSeq)]]
        bestTheta = glist
      }


      tmp = glist$a
      tmp[tmp == 0] = 1e300
      whichMin = which.min(tmp)
      glist$a[whichMin] = 0
      cat("knz =", knz, " --- ")
      # knz = knz - 1L
      knz = sum(glist$a > 0)
      # cat("knz =", knz, " --- ")
      glist$a = glist$a / sum(glist$a)
      # denMat = computeDensityMat(X, glist, withAlpha = F)
    }


    tmp = which(bestTheta$a > 0)
    bestTheta$a = bestTheta$a[tmp]
    bestTheta$mu = bestTheta$mu[, tmp]
    bestTheta$sigma = bestTheta$sigma[, tmp]


    denWithAlpha = computeDensityMat(X, bestTheta, withAlpha = T)
    tmp = apply(denWithAlpha, 1, function(x) which.max(x))
    clusterMember = aggregate(list(1L : nrow(X)), list(tmp), function(x) x)[[2]]
    for(i in 1L : length(clusterMember)) clusterMember[[i]] = sort(clusterMember[[i]])


    list(bestTheta = bestTheta, fitted = apply(computeDensityMat(X, bestTheta, withAlpha = T), 1, function(x) sum(x)), minLoss = minLoss, clusterMember = clusterMember)
  }


}




for(k in 1L : 1000){ cat(k, "")


g = 25
# dat = t(t(iris[1:4]) / apply(iris[1:4], 2, function(x) sd(x)))
dat = as.matrix(iris[1:4])


dimnames(dat) = NULL
mu = t(dat[sample(1L : nrow(iris), g), ])
dimnames(mu) = NULL
a = rep(1 / g, g)
r = 0.1
s = matrix(rep(as.numeric(diag(max(diag(var(dat))) * r, nrow = 4, ncol = 4)), g), ncol = g)


glist = list(a = a, mu = mu, sigma = s)
# sink("rst.txt")
# rst = FJgmmOriginal(dat, glist, kmin = 2L, eps = 1e-5)
# sink()
# sink("rst2.txt")
# rst2 = FJgmm(dat, glist, kmin = 2L, eps = 1e-5)
# str(rst2)
# sink()
sink("rst3.txt")
rst3 = FJgmmMatlab(dat, glist, kmin = 2L, eps = 1e-5)
sink()
str(rst3)


sink("debug.txt")
rst4 = GMKMcharlie::GMfj(t(dat), alpha = glist$a, mu = glist$mu, sigma = glist$sigma, maxCore = 1, tlimit = 100)
sink()
str(rst3)
str(rst4)
# str(rst)
# str(rst2)
# str(rst3)


tmp = range(c(range(unlist(rst3) / unlist(rst4) - 1)))
# tmp = range(c(range(unlist(rst) / unlist(rst2) - 1), range(unlist(rst) / unlist(rst3) - 1)))
str(tmp)


if(max(abs(tmp)) > 1e-5) break
}
save.image()
write("", file = paste0("kEquals", k))
print("k = ")
print(k)




# Check against FJ matlab code
if(T)
{


  dat = as.matrix(read.csv("C:/Users/i56087/Downloads/mixturecode2/data.csv", header = F))
  dimnames(dat) = NULL
  N = ncol(dat)
  d = nrow(dat)
  kmax = 25
  glist = list()
  glist$a = unlist(read.csv("C:/Users/i56087/Downloads/mixturecode2/initialestpp.csv", header = F))
  names(glist$a) = NULL
  glist$mu = as.matrix(read.csv("C:/Users/i56087/Downloads/mixturecode2/initialmu.csv", header = F))
  dimnames(glist$mu) = NULL
  glist$sigma = matrix(rep(as.numeric(as.matrix(read.csv("C:/Users/i56087/Downloads/mixturecode2/cov.csv", header = F))), kmax), ncol = kmax)
  rst = FJgmmMatlab(t(dat), glist, kmin = 2L, eps = 1e-5)
  str(rst)


}

























# a = glist$a
# m = glist$mu
# s = glist$sigma



while(T)
{
  initialG = 20L
  N = 10000L
  d = 10L
  underlyingG = 4L
  dat = matrix(unlist(lapply(1L : underlyingG, function(x)
  {
    matrix(rnorm(as.integer(N / underlyingG) * d), nrow = d) * runif(1, 0.9, 1.1) + runif(d, -10, 10)
  })), nrow = d)
  N = ncol(dat)


  mu = dat[, sample(ncol(dat), initialG)]
  s = matrix(rep(var(t(dat)), initialG), nrow = d * d) / initialG
  a = rep(1 / initialG, initialG)


  # sink("debug.txt")
  # input data matrix should be column-wise, while mu and sigma's numbers of rows equal d and d ^ 2.
  system.time({fittedGMr = FJgmm(t(dat), list(a = a, mu = mu, sigma = s), kmin = 2L, eps = 1e-5)})
  # sink()


  # sink("debug2.txt")
  system.time({fittedGMcpp = GMKMcharlie::paraGmmFJ(
    dat, Xw = numeric(0), Gmin = 2, alpha = a, mu = mu, sigma = s, maxCore = 7,
    maxIter = 1000, eigenRatioLim = 0, convergenceEPS = 1e-5, alphaEPS = 1e-16,
    tlimit = 3600, verbose = 0)})
  # sink()
  cat(abs(fittedGMr$minLoss / fittedGMcpp$loss - 1), "")


  if(abs(fittedGMr$minLoss / fittedGMcpp$loss - 1) > 1e-2) break
}




g = 4L
dat = t(t(iris[1:4]) / apply(iris[1:4], 2, function(x) sd(x)))
dimnames(dat) = NULL
mu = t(dat[sample(1L : nrow(iris), g), ])
# mu = t(iris[sample(1L : nrow(iris), g), 1:4])
# mu = GMKMcharlie:::findSpreadedMeanWrapper(t(iris[1:4]), g, maxCore = 7)
dimnames(mu) = NULL
# a = runif(g); a = a / sum(a)
a = rep(1 / g, g)
# s = matrix(rep(as.numeric(var(iris[1:4])), g), ncol= g)
r = 1
s = matrix(rep(as.numeric(diag(mean(diag(var(dat))) * r, nrow = 4, ncol = 4)), g), ncol = g)





# Test density calculation. There is no problem in density calculations
if(F)
{
  denMat = computeDensityMat(as.matrix(iris[1:4]), list(a = a, mu = mu, sigma = s), withAlpha = F)
  dimnames(denMat) = NULL
  tmp = as.matrix(as.data.frame(mapply(function(m, s)
  {
    apply(t(iris[1:4]), 2, function(x)
    {
      GMKMcharlie::testGdensity(x, mu = m, sigma = s)
    })
  }, as.data.frame(mu), as.data.frame(s), SIMPLIFY = F)))
  dimnames(tmp) = NULL
}





























if(F)
{
  sink("debug3.txt")
  fittedGM = GMKMcharlie::paraGmmFJ(
    t(iris[1:4]) / apply(iris[1:4], 2, function(x) sd(x)), G = 20, Gmin = 2, maxCore = 7,
    maxIter = 1000, eigenRatioLim = 0, convergenceEPS = 1e-5, alphaEPS = 1e-16,
    tlimit = 3600, verbose = 10)
  sink()
}



length(tmp$bestTheta$a); length(fittedGM$alpha)
tmp$minLoss; fittedGM$loss


range(fittedGM$fitted / tmp$fitted - 1)
# sort(fittedGM$alpha); sort(tmp$bestTheta$a)
odr = order(tmp$fitted)
plot(tmp$fitted[odr], col = "darkblue")
lines(fittedGM$fitted[odr], lwd = 0.1, col = "red")












