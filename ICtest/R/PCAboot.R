##################################################################################
# Bootstrapping PCA
##################################################################################

# internal
# Function to multiply a vector with a random orthogonal matrix

prodRandO <- function(x, k) x %*% rorth(k)

# internal function to compute test statistic

PCA_boot_teststatistic <- function(X, k, S1, Sargs)
    {
    n <- nrow(X)
    p <- ncol(X)

    LocScat <- do.call(S1, c(list(X), Sargs))
    MEAN <- LocScat[[1]]
    COV <- LocScat[[2]]

    Xc <- sweep(X, 2, MEAN, "-")
    EV.COV <- eigen(COV, symmetric=TRUE)

    W <- EV.COV$vectors
    D <- EV.COV$values

    ind <- (k+1):p
    EVpk <- D[ind]
    meanEVpk <- mean(EVpk)

    TEST.STATISTIC.X <-  sum((EVpk - meanEVpk)^2)/(meanEVpk^2*(p-k))
    return(TEST.STATISTIC.X)
    }

# Internal
# Function to do the bootstrap sample and compute the corresponding test statistic

PCA_boot_orth_s1 <- function(Z1, Z2, S1, W, MEAN, k, Sargs)
    {
    n <- nrow(Z2)
    p2 <- ncol(Z2)

    ind <- sample(1:n, replace=TRUE)

    Z2tilde <- t(apply(Z2[ind,], 1, prodRandO, k=p2))
    Ztilde <- cbind(Z1[ind,],Z2tilde)
    Zstar <- Ztilde

    Xstar <- tcrossprod(Zstar, W) #Zstar %*% t(W)

    Xstar <- sweep(Xstar, 2, MEAN, "+")

    TEST.STATISTIC.Xstar <-  PCA_boot_teststatistic(X=Xstar, k=k, S1=S1, Sargs=Sargs)

    TEST.STATISTIC.Xstar
    }

# internal
# General function for PCA bootstrapping using strategy 1

PCA_subsphere_boot_s1 <- function(X, k, S1, n.boot, Sargs)
    {
    n <- nrow(X)
    p <- ncol(X)
    LocScat <- do.call(S1, c(list(X), Sargs))
    MEAN <- LocScat[[1]]
    COV <- LocScat[[2]]
    Xc <- sweep(X, 2, MEAN, "-")
    EV.COV <- eigen(COV, symmetric=TRUE)
    W <- EV.COV$vectors
    D <- EV.COV$values
    ind <- (k+1):p
    EVpk <- D[ind]
    meanEVpk <- mean(EVpk)
    TEST.STATISTIC.X <-  sum((EVpk - meanEVpk)^2)/(meanEVpk^2*(p-k))
    names(TEST.STATISTIC.X) = "T"
    PARAMETER <- n.boot
    names(PARAMETER) <- c("replications")
    Z <- Xc %*% W
    # Column names for the variables
    colnames(Z) <- paste0("PC.",1:p)
    Z1 <- Z[ ,1:k, drop=FALSE]
    Z2 <- Z[ ,ind, drop=FALSE]

    TEST.STATISTICS.Xstar <- t(replicate(n.boot, PCA_boot_orth_s1(Z1, Z2, S1, W, MEAN, k, Sargs=Sargs)))
    
    PVAL <- (sum(TEST.STATISTIC.X<=TEST.STATISTICS.Xstar)+1)/(n.boot+1)
    ALTERNATIVE <- paste0("the last ",p-k, " eigenvalues are not equal")
    RES <- list(statistic = n*TEST.STATISTIC.X, p.value=PVAL, parameter=PARAMETER, alternative = ALTERNATIVE, k=k, W=W, S=Z, D=D, MU=MEAN, SCATTER=COV)
    RES
    }

#########################################################################################

# Strategy 2

# Function to do the bootstrap sample and compute the corresponding test statistic

PCA_boot_orth_s2 <- function(Z, S1, W, LAMBDA, MEAN, k, Sargs)
    {
    n <- nrow(Z)
    p <- ncol(Z)
    ind <- sample(1:n, replace=TRUE)
    Ztilde <- t(apply(Z[ind,], 1, prodRandO, k=p))
    Xstar <- sweep(Ztilde,2,LAMBDA,"*")
    Xstar <- Xstar %*% t(W)
    Xstar <- sweep(Xstar, 2, MEAN, "+")
    TEST.STATISTIC.Xstar <-  PCA_boot_teststatistic(X=Xstar, k=k, S1=S1, Sargs=Sargs)
    c(TEST.STATISTIC.Xstar)
    }


# Internal
# General function for PCA bootstrapping using strategy 2

PCA_subsphere_boot_s2 <- function(X, S1,k, n.boot, Sargs)
    {
    n <- nrow(X)
    p <- ncol(X)
    LocScat <- do.call(S1, c(list(X), Sargs))
    MEAN <- LocScat[[1]]
    COV <- LocScat[[2]]
    Xc <- sweep(X, 2, MEAN, "-")
    EV.COV <- eigen(COV, symmetric=TRUE)
    W <- EV.COV$vectors
    D <- EV.COV$values
    ind <- (k+1):p
    EVpk <- D[ind]
    meanEVpk <- mean(EVpk)
    TEST.STATISTIC.X <-  sum((EVpk - meanEVpk)^2)/(meanEVpk^2*(p-k))
    names(TEST.STATISTIC.X) = "T"
    PARAMETER <- n.boot
    names(PARAMETER) <- c("replications")
    LAMBDA <- sqrt(c(D[1:k], rep(meanEVpk,p-k)))
    #LAMBDA.inv <- 1/LAMBDA
    LAMBDA.inv <- 1/D
    Z <- Xc %*% W
    S <- Z
    colnames(S) <- paste0("PC.",1:p)
    Z <- sweep(Z,2,LAMBDA.inv,"*")
    TEST.STATISTICS.Xstar <- t(replicate(n.boot, PCA_boot_orth_s2(Z, S1, W, LAMBDA, MEAN, k, Sargs=Sargs)))
    PVAL <- (sum(TEST.STATISTIC.X<=TEST.STATISTICS.Xstar)+1)/(n.boot+1)
    ALTERNATIVE <- paste0("the last ",p-k, " eigenvalues are not equal")
    RES <- list(statistic = n*TEST.STATISTIC.X, p.value=PVAL, parameter=PARAMETER, alternative = ALTERNATIVE, k=k, W=W, S=S, D=D, MU=MEAN, SCATTER=COV)
    RES
    }

PCAboot <- function(X, k, n.boot = 200, s.boot = "B1" , S = MeanCov, Sargs=NULL)
    {

    DNAME <- deparse(substitute(X))
    SNAME <- deparse(substitute(S))
    
    X <- as.matrix(X)

    s.boot <- match.arg(s.boot, c("B1", "B2"))
    if (s.boot == "B1") {
       RES <- PCA_subsphere_boot_s1(X=X, S1=S, k=k, n.boot=n.boot, Sargs=Sargs)
       METHOD <- paste0("PCA subsphericity bootstrapping test using ", SNAME, " and strategy B1")
       } else {
       RES <- PCA_subsphere_boot_s2(X=X, S1=S, k=k, n.boot=n.boot, Sargs=Sargs)
       METHOD <- paste0("PCA subsphericity bootstrapping test using ", SNAME, " and strategy B2")
       }

    RES <- c(RES, method=METHOD, data.name=DNAME, s.boot=s.boot, scatter=SNAME)
    # Reorder the list
    RES <- RES[c(1:3, 11:12, 4:10, 14, 13)]
    class(RES) <- c("ictest", "htest")
    RES
    }
