PCAschott <- function(X, k)
    {
    DNAME <- deparse(substitute(X))
    n <- nrow(X)
    p <- ncol(X)
    X <- as.matrix(X)

    MEAN <- colMeans(X)
    Xc <- sweep(X,2,MEAN,"-")
    COV <- crossprod(Xc)/(n-1)


    METHOD <- "PCA subspericity test using Schott's test"

    EV <- eigen(COV, symmetric = TRUE)
    EV.VALUES <- EV$values
    ind <- (k+1):p
    MEANev <- mean(EV.VALUES[ind])

    cn <- (n-1-k) - 1/6*(2*(p-k) + 1 + (p-k)/2) + MEANev^2* 1/sum((EV.VALUES[-ind]-MEANev)^2)

    TESTSTATISTIC <- cn * ((p-k)*log(MEANev)- sum(log(EV.VALUES[ind])))
    names(TESTSTATISTIC) = "T"
    PARAMETER <- ((p-k)*(p-k+1)/2)-1
    names(PARAMETER) <- c("df")
    PVAL <- 1- pchisq(TESTSTATISTIC, df=PARAMETER)
    ALTERNATIVE <- paste0("the last ",p-k, " eigenvalues are not equal")

    W <- EV$vectors
    S <- tcrossprod(Xc,W)
    colnames(S) <- paste0("PC.",1:p)

    RES <- list(statistic = TESTSTATISTIC, p.value = PVAL, parameter = PARAMETER, method=METHOD, data.name=DNAME, alternative = ALTERNATIVE, k=k, W=W, S=S, D=EV.VALUES,
    MU=MEAN, SCATTER=COV)
    class(RES) <- c("ictest", "htest")
    RES
    }
