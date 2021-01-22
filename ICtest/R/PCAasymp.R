PCAasymp <- function(X, k, scatter="cov", ...)
    {
    DNAME <- deparse(substitute(X))
    n <- nrow(X)
    p <- ncol(X)
    X <- as.matrix(X)

    scatter <-  match.arg(scatter, c("cov", "tyler"))

    if (scatter=="cov"){
            MEAN <- colMeans(X)
            Xc <- sweep(X,2,MEAN,"-")
            COV <- crossprod(Xc)/n
            MAHA <- rowSums(Xc %*% solve(COV) * Xc)
            sigma1 <- (mean(MAHA^2))/(p*(p+2))
            METHOD <- "PCA subspericity test using cov"
            } else {
            TYLER <- HR.Mest(X, ...)
            MEAN <- TYLER$center
            COV <-  TYLER$scatter
            Xc <-  sweep(X,2,MEAN,"-")
            sigma1 <- (p+2)/p
            METHOD <- "PCA subspericity test using Tyler's shape matrix"
            }

    EV <- eigen(COV, symmetric = TRUE)
    EV.VALUES <- EV$values
    ind <- (k+1):p
    MEANev <- mean(EV.VALUES[ind])
    VAR.K <- sum((EV.VALUES[ind]-MEANev)^2)/(p-k)
    TESTSTATISTIC <- n*VAR.K*(p-k)/(2*sigma1*MEANev^2)
    names(TESTSTATISTIC) = "T"
    PARAMETER <- (p-k-1)*(p-k+2)/2
    names(PARAMETER) <- c("df")
    PVAL <- 1- pchisq(TESTSTATISTIC, df=PARAMETER)
    ALTERNATIVE <- paste0("the last ",p-k, " eigenvalues are not equal")

    W <- EV$vectors
    # S <- tcrossprod(Xc,W)
    S <- Xc %*% W
    colnames(S) <- paste0("PC.",1:p)

    RES <- list(statistic = TESTSTATISTIC, p.value = PVAL, parameter = PARAMETER, method=METHOD, data.name=DNAME, alternative = ALTERNATIVE, k=k, W=W, S=S, D=EV.VALUES, 
    MU=MEAN, SCATTER=COV, sigma1 = sigma1)
    class(RES) <- c("ictest", "htest")
    RES
    }
