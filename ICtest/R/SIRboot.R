# internal:
# Function to compute the test statistic

SIR_boot_teststatistic <- function(X, y, k, h, ...)
    {
    n <- nrow(X)
    p <- ncol(X)
    MEAN <- colMeans(X)
    Xc <- sweep(X, 2, MEAN, "-")
    COV <- crossprod(Xc)/n

    EVD.COV <- eigen(COV, symmetric=TRUE)
    COV.inv.sqrt <- EVD.COV$vectors %*% diag((1/EVD.COV$values)^0.5) %*%  t(EVD.COV$vectors)

    Y <- tcrossprod(Xc, COV.inv.sqrt)
    S2 <- covSIR(Y, y, h=h, ...)
    EVD.S2 <- eigen(S2, symmetric = TRUE)

    D <- EVD.S2$values

    W <- crossprod(EVD.S2$vectors, COV.inv.sqrt)

    TEST.STATISTIC.X <- sum(D[(k+1):p])
    return(TEST.STATISTIC.X)
    }

# internal:
# Function to create the bootstrap sample and to compute its test statistic

SIR_boot_normal_all <- function(Z1,Z2,y, Winv, MEAN, k, h, ...)
    {
    n <- nrow(Z1)

    ind <- sample(1:n, replace=TRUE)

    Z1tilde <- Z1[ind,]
    Z2tilde<-Z2[sample(1:n,replace=TRUE),]
    Zstar <- cbind(Z1tilde,Z2tilde)

    Xstar <- tcrossprod(Zstar, Winv)
    Xstar <- sweep(Xstar, 2, MEAN, "+")

    ystar <- y[ind]
    TEST.STATISTIC.Xstar <- SIR_boot_teststatistic(Xstar, ystar, k, h=h, ...)
    TEST.STATISTIC.Xstar
    }

# General function for SIR bootstrapping:

SIRboot <- function(X, y, k, h=10, n.boot = 200, ...)
    {
    DNAME <- deparse(substitute(X))
    n <- nrow(X)
    p <- ncol(X)
    MEAN <- colMeans(X)
    Xc <- sweep(X, 2, MEAN, "-")
    COV <- crossprod(Xc)/n

    EVD.COV <- eigen(COV, symmetric=TRUE)
    COV.inv.sqrt <- EVD.COV$vectors %*% diag((1/EVD.COV$values)^0.5) %*%  t(EVD.COV$vectors)

    Y <- tcrossprod(Xc, COV.inv.sqrt)
    S2 <- covSIR(Y, y, h=h,...)
    EVD.S2 <- eigen(S2, symmetric = TRUE)

    D <- EVD.S2$values
    W <- crossprod(EVD.S2$vectors, COV.inv.sqrt)

    TEST.STATISTIC.X <-  sum(D[(k+1):p])

    Z <- tcrossprod(Xc, W)

    Z1 <- Z[,0:k, drop=FALSE]
    Z2 <- Z[,(k+1):p, drop=FALSE]

    Winv <- solve(W)

    TEST.STATISTICS.Xstar <- replicate(n.boot, SIR_boot_normal_all(Z1, Z2, y,  Winv, MEAN, k, h, ...))
    PVAL <- (sum(TEST.STATISTIC.X<TEST.STATISTICS.Xstar)+1)/(n.boot+1)
    METHOD <- "SIR bootstrapping test for subspace dimension"
    ALTERNATIVE <- paste0("the last ",p-k, " eigenvalues are not zero")
    PARAMETER <- n.boot
    names(PARAMETER) <-  "replications"
    names(TEST.STATISTIC.X) <- "T"
    
    colnames(Z) <- paste0("SIC.",1:p)
    
    RES <- list(statistic = TEST.STATISTIC.X, p.value = PVAL, parameter = PARAMETER, method=METHOD, data.name=DNAME, alternative = ALTERNATIVE, k=k, W=W, S=Z, D=D,
    MU=MEAN)
    class(RES) <- c("ictest", "htest")
    RES
    }
