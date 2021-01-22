
# ICA test for subnormality

# internal
# Function to make a bootstrap sample for each column of a matrix

SampFunction <- function(x,n)  x[sample(1:n, replace=TRUE)]

# internal
# Function to compute ICA test statistic for a sample X and hypothesis k.

ICA_boot_teststatistic <- function(X, k)
    {
    n <- nrow(X)
    p <- ncol(X)
    MEAN <- colMeans(X)
    Xc <- sweep(X, 2, MEAN, "-")
    COV <- crossprod(Xc)/n

    EVD.COV <- eigen(COV, symmetric=TRUE)
    COV.inv.sqrt <- EVD.COV$vectors %*% diag((1/EVD.COV$values)^0.5) %*%  t(EVD.COV$vectors)

    Y <- tcrossprod(Xc, COV.inv.sqrt)
    r <- sqrt(rowSums(Y^2))
    Y <- r * Y
    COV4 <- crossprod(Y)/n
    EVD.COV4 <- eigen(COV4, symmetric = TRUE)

    D <- EVD.COV4$values
    orderD <- order((D-(p+2))^2, decreasing = TRUE)
    Dord <- D[orderD]
    #W <- crossprod(EVD.COV4$vectors, COV.inv.sqrt)[orderD,]

    TEST.STATISTIC.X <-  sum((Dord[(k+1):p]-(p+2))^2)
    return(TEST.STATISTIC.X)
    }

# Internal
# Function to create the bootstrap sample and compute the test statistic using Statregy 1

ICA_boot_normal_S1 <- function(Z1, Winv, MEAN, k)
    {
    n <- nrow(Z1)
    p <- ncol(Winv)

    Z2tilde <- matrix(rnorm(n*(p-k)), nrow=n)
    Z1tilde <- Z1[sample(1:n, n, replace=TRUE), , drop=FALSE]
    Zstar <- cbind(Z1tilde,Z2tilde)

    Xstar <- tcrossprod(Zstar, Winv)
    Xstar <- sweep(Xstar, 2, MEAN, "+")
    TEST.STATISTIC.Xstar <- ICA_boot_teststatistic(Xstar, k)

    TEST.STATISTIC.Xstar
    }

# interal
# Function to create the bootstrap sample and compute the test statistic using Statregy 2

ICA_boot_normal_S2 <- function(Z1, Winv, MEAN, k)
    {
    n <- nrow(Z1)
    p <- ncol(Winv)

    Z2tilde <- matrix(rnorm(n*(p-k)), nrow=n)
    Z1tilde <- apply(Z1,2,SampFunction, n=n)
    Zstar <- cbind(Z1tilde,Z2tilde)

    Xstar <- Xstar2 <- tcrossprod(Zstar, Winv)
    Xstar <- sweep(Xstar, 2, MEAN, "+")
    TEST.STATISTIC.Xstar <- ICA_boot_teststatistic(Xstar, k)

    TEST.STATISTIC.Xstar
    }

# General ICA function using Strategy 1

ICA_subnormal_boot_S1 <- function(X,k, n.boot = 200)
    {
    n <- nrow(X)
    p <- ncol(X)
    MEAN <- colMeans(X)
    Xc <- sweep(X, 2, MEAN, "-")
    COV <- crossprod(Xc)/n

    EVD.COV <- eigen(COV, symmetric=TRUE)
    COV.inv.sqrt <- EVD.COV$vectors %*% diag((1/EVD.COV$values)^0.5) %*%  t(EVD.COV$vectors)

    Y <- tcrossprod(Xc, COV.inv.sqrt)
    r <- sqrt(rowSums(Y^2))
    Y <- r * Y
    COV4 <- crossprod(Y)/n
    EVD.COV4 <- eigen(COV4, symmetric = TRUE)

    D <- EVD.COV4$values
    orderD <- order((D-(p+2))^2, decreasing = TRUE)
    Dord <- D[orderD]
    W <- crossprod(EVD.COV4$vectors, COV.inv.sqrt)[orderD,]


    TEST.STATISTIC.X <-  sum((Dord[(k+1):p]-(p+2))^2)
    names(TEST.STATISTIC.X) = "T"
    PARAMETER <- n.boot
    names(PARAMETER) <- c("replications")

    Z <- tcrossprod(Xc, W)
    
    Z1 <- Z[,0:k, drop=FALSE]

    Winv <- solve(W)

    TEST.STATISTICS.Xstar <- t(replicate(n.boot, ICA_boot_normal_S1(Z1,  Winv, MEAN, k)))
    
    names(Z) <- paste0("IC.",1:p)
    PVAL <- (sum(TEST.STATISTIC.X<TEST.STATISTICS.Xstar)+1)/(n.boot+1)
    ALTERNATIVE <- paste0("the last ",p-k, " components are not gaussian")
    RES <- list(statistic = n*TEST.STATISTIC.X, p.value=PVAL, parameter=PARAMETER, alternative = ALTERNATIVE, k=k, W=W, S=Z, D=D, MU=MEAN)
    RES
    }



# General ICA function using Strategy 2

ICA_subnormal_boot_S2 <- function(X,k, n.boot = 200)
    {
    n <- nrow(X)
    p <- ncol(X)
    MEAN <- colMeans(X)
    Xc <- sweep(X, 2, MEAN, "-")
    COV <- crossprod(Xc)/n

    EVD.COV <- eigen(COV, symmetric=TRUE)
    COV.inv.sqrt <- EVD.COV$vectors %*% diag((1/EVD.COV$values)^0.5) %*%  t(EVD.COV$vectors)

    Y <- tcrossprod(Xc, COV.inv.sqrt)
    r <- sqrt(rowSums(Y^2))
    Y <- r * Y
    COV4 <- crossprod(Y)/n
    EVD.COV4 <- eigen(COV4, symmetric = TRUE)

    D <- EVD.COV4$values
    orderD <- order((D-(p+2))^2, decreasing = TRUE)
    Dord <- D[orderD]
    W <- crossprod(EVD.COV4$vectors, COV.inv.sqrt)[orderD,]


    TEST.STATISTIC.X <-  sum((Dord[(k+1):p]-(p+2))^2)
    names(TEST.STATISTIC.X) = "T"
    PARAMETER <- n.boot
    names(PARAMETER) <- c("replications")

    Z <- tcrossprod(Xc, W)
    colnames(Z) <- paste0("IC.",1:p)
    
    Z1 <- Z[,0:k, drop=FALSE]

    Winv <- solve(W)

    TEST.STATISTICS.Xstar <- t(replicate(n.boot, ICA_boot_normal_S2(Z1,  Winv, MEAN, k)))

    PVAL <- (sum(TEST.STATISTIC.X<TEST.STATISTICS.Xstar)+1)/(n.boot+1)
    ALTERNATIVE <- paste0("the last ",p-k, " components are not gaussian")
    RES <- list(statistic = n*TEST.STATISTIC.X, p.value=PVAL, parameter=PARAMETER, alternative = ALTERNATIVE, k=k, W=W, S=Z, D=D, MU=MEAN)
    RES
    }


FOBIboot <- function(X, k, n.boot = 200, s.boot = "B1")
    {

    DNAME <- deparse(substitute(X))
    X <- as.matrix(X)

    s.boot <- match.arg(s.boot, c("B1", "B2"))
    if (s.boot == "B1") {
       RES <- ICA_subnormal_boot_S1(X=X,k=k,n.boot)
       METHOD <- "ICA subgaussianity bootstrapping test using  FOBI and strategy B1"
       } else {
       RES <- ICA_subnormal_boot_S2(X,k,n.boot)
       METHOD <- "ICA subgaussianity bootstrapping test using  FOBI and strategy B2"
       }

    RES <- c(RES, method=METHOD, data.name=DNAME, s.boot=s.boot)
    # Reorder the return list
    RES <- RES[c(1:3, 10:11, 4:9, 12)]
    class(RES) <- c("ictest", "htest")
    RES
    }
