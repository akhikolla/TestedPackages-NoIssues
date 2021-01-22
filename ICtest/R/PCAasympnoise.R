PCAasympnoise <- 
  function (X, k, noise=c("Gaussian","Elliptical"))
  {
    DNAME <- deparse(substitute(X))
    noise <- match.arg(noise)
    n <- nrow(X)
    p <- ncol(X)
    X <- as.matrix(X)
    
    MEAN <- colMeans(X)
    Xc <- sweep(X, 2, MEAN, "-")
    COV <- crossprod(Xc)/n
    METHOD <- "PCA dimension test in a noisy model"
    
    EV <- eigen(COV, symmetric = TRUE)
    EV.VALUES <- EV$values
    ind <- (k + 1):p
    MEANev <- mean(EV.VALUES[ind])
    VAR.K <- sum((EV.VALUES[ind] - MEANev)^2)/(p - k)
    
    W <- EV$vectors
    S <- Xc %*% W
    
    if (noise=="Gaussian") {sigma1 <- 1} else {
      S2 <- S[, ind]
      MAHA <- rowSums(S2 %*% solve(cov(S2)) * S2)
      sigma1 <- (mean(MAHA^2))/((p-k) * (p-k + 2))
    }
    
    
    TESTSTATISTIC <- n * VAR.K * (p - k)/(2 * sigma1 * MEANev^2)
    names(TESTSTATISTIC) = "T"
    PARAMETER <- (p - k - 1) * (p - k + 2)/2
    names(PARAMETER) <- c("df")
    PVAL <- 1 - pchisq(TESTSTATISTIC, df = PARAMETER)
    ALTERNATIVE <- paste0("the last ", p - k, " eigenvalues are not equal")
    W <- EV$vectors
    S <- Xc %*% W
    colnames(S) <- paste0("PC.", 1:p)
    RES <- list(statistic = TESTSTATISTIC, p.value = PVAL, parameter = PARAMETER, 
                method = METHOD, data.name = DNAME, alternative = ALTERNATIVE, 
                k = k, W = W, S = S, D = EV.VALUES, MU = MEAN, SCATTER = COV, 
                sigma1 = sigma1)
    class(RES) <- c("ictest", "htest")
    RES
  }