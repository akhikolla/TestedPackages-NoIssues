mcfa<- function(Y, g, q, itmax=50, nkmeans=20, nrandom=20, tol=1.e-5,
                initClust=NULL, initMethod='eigenA', convMeas='diff',
                errorMsg=FALSE, verbose=FALSE, ...){
    
    if(!is.matrix(Y)) Y <- as.matrix(Y)
    p <- ncol(Y)
    if(is.null(p))
        stop("the data must have more than one variable")
    #return(ERRMSG <- "Error: the data must have more than one variable")
    n <- nrow(Y)
    if(p < q)
        stop("the number of factors must not be greater
    than the number of variables")

    startClust <- Starts(Y, g, initClust, nkmeans, nrandom)
    maxinit <- ncol(startClust)
    if(is.null(maxinit)) maxinit <- 1
    
    ERRMSG <- NULL
    maxLOGL  <- -Inf
    for (ii in seq_len(maxinit)) {
        if(min(table(startClust[, ii]) == 1)){
            when <- paste("At start", ii)
            what <- "Initial partition was not used as it has a cluster of one
            sample."
            ERRMSG <- rbind(ERRMSG, cbind(when, what))
            next
        }
        startModel <- try(init.est.para.mcfa(Y, g, q, startClust[, ii],
            initMethod))
        if(class(startModel)[1] == "try-error"){
            when <- paste("At start", ii)
            what <- "Failed to estimate initial parameters"
            ERRMSG <- rbind(ERRMSG, cbind(when, what))
            next
        }
        estModel <- est.mcfa(Y, g, q, itmax, tol, startModel$pivec, 
            startModel$A,startModel$xi, startModel$omega, 
            startModel$D, convMeas)
        if((class(estModel)[1] == "mcfa")){
            if (estModel$logL > maxLOGL){
                Hmodel <- estModel
                maxLOGL <- Hmodel$logL
            }
            if(verbose){
                message(sprintf("g = %i, q = %i, initialization %i logL %8.4f,
                    maxlogL = %8.4f \n",
                    g, q, ii, estModel$logL, maxLOGL))
            }
        }
        if(class(estModel)[1] == "error"){
            when <- paste("At start", ii)
            what <- estModel
            ERRMSG <- rbind(ERRMSG, cbind(when, what))
        }
    }
    
    if(!exists("Hmodel")){
        stop("Error: Failed to Estimate a Model. See Error Messages.")
        return(ERRMSG)
    }
    
    CH <- chol(t(Hmodel$A) %*% Hmodel$A)
    Hmodel$A <- Hmodel$A %*% solve(CH)
    Hmodel$xi <- CH %*% Hmodel$xi
    for(i in seq_len(g)){
        Hmodel$omega[,,i] <- CH %*% Hmodel$omega[,,i] %*% t(CH)
    }
    d <- (g-1) + p + q*(p+g) + g*q*(q+1)/2 - q*q
    Hmodel$tau <- do.call('tau.mcfa', c(list(Y=Y), Hmodel))
    Hmodel$BIC <- -2*Hmodel$logL + d*log(n)
    Hmodel$clust <- apply(Hmodel$tau, 1, which.max)
    Hmodel <- append(Hmodel, do.call('factor.scores.mcfa',
        c(list(Y=Y), Hmodel)))
    Hmodel$call <- match.call()
    if(errorMsg == TRUE)
        Hmodel$ERRMSG <- ERRMSG
    class(Hmodel) <- "mcfa"
    return(Hmodel)
}

est.mcfa <- function(Y, g, q, itmax, tol, pivec, 
        A, xi, omega, D, convMeas, ...){
    if(!is.matrix(Y)) Y <- as.matrix(Y)
    p <- ncol(Y); if(is.null(p)) p <- 1
    n <- nrow(Y)
    
    fit <- list(g=g, q=q, pivec=pivec, A=A, xi=xi, omega=omega, D=D)
    
    fit$logL <- loglike.mcfa(Y, g, q, pivec, A, xi, omega, D)
    if(class(fit$logL)[1] == 'character'){
        FIT <- paste('Failed in computing the 
            log-likelihood before the EM-steps,', fit$logL)
        class(FIT) <- "error"
        return(FIT)
    }
    for(niter in seq_len(itmax)) {
        FIT <- do.call('Mstep.mcfa', c(list(Y=Y), fit))
        if(class(FIT)[1] == 'error'){
            FIT <- paste('Computational error in ', niter,
                'iteration of the M-step:',FIT)
            class(FIT) <- "error"
            return(FIT)
        }
        FIT$logL <- do.call('loglike.mcfa', c(list(Y=Y), FIT))
        if(class(FIT$logL)[1] == 'character'){
            FIT <- paste('Failed to compute log-likelihood after ', niter,
                'th the M-step:', FIT$logL, sep='')
            class(FIT) <- "error"
            return(FIT)
        }
        if ((FIT$logL == -Inf) | is.na(FIT$logL)){
            FIT <- paste('Log likelihood computed after the', niter,
                'th iteration of the M-step is not finite', sep='')
            class(FIT) <- "error"
            return(FIT)
        }
        if((convMeas == "diff") & (abs(FIT$logL-fit$logL) < tol))
            break
        if((convMeas == "ratio") & (abs((FIT$logL-fit$logL)/FIT$logL) < tol))
            break
        fit <- FIT
    }
    class(FIT) <- "mcfa"
    return(FIT)
}

loglike.mcfa <- function(Y, g, q, pivec, A, xi, omega, D, ...){
    p <- ncol(Y); if(is.null(p)) p <- 1
    n <- nrow(Y)
    
    Fji <- array(NA, c(n,g))
    InvD <- diag(1/diag(D))
    for(i in seq_len(g)){
        InvS <- try(InvD - InvD %*% A %*%
            chol.inv(chol.inv(omega[,,i]) + t(A) %*% InvD %*% A) %*%
            t(A) %*% InvD)
        if(class(InvS)[1] == "try-error") {
            return(loglike <- paste('ill-cond. or sing. Sigma_', i,sep=''))
        }
        logdetD <- log(det(as.matrix(omega[,,i]))) + sum(log(diag(D))) +
            log(det(chol.inv(omega[,,i]) + t(A) %*% InvD %*% A))
        
        
        
        MhalDist <- stats::mahalanobis(Y, t(A %*% xi[, i, drop=FALSE]), InvS, 
            inverted=TRUE)
        Fji[,i] <- -0.5*MhalDist - (p/2)*log(2*pi) - 0.5*logdetD
    }
    Fji <- sweep(Fji, 2, log(pivec), '+')
    Fjmax <- apply(Fji, 1, max)
    Fji <- sweep(Fji, 1, Fjmax, '-')
    loglike <- sum(Fjmax, log(rowSums(exp(Fji))))
    return(loglike)
}

tau.mcfa <- function(Y, g, q, pivec, A, xi, omega, D, ...){
    p <- ncol(Y); if(is.null(p)) p <- 1
    n <- nrow(Y)
    Fji <- array(NA, c(n,g))
    InvD <- diag(1/diag(D))
    for(i in seq_len(g)){
        InvS <- InvD - InvD %*% A %*%
            chol.inv(chol.inv(omega[,,i]) + t(A) %*% InvD %*% A) %*%
            t(A) %*% InvD
        logdetD <- log(det(as.matrix(omega[,,i]))) + sum(log(diag(D))) +
            log(det(chol.inv(omega[,,i]) + t(A) %*% InvD %*% A))
        MhalDist <- stats::mahalanobis(Y, t(A %*% xi[, i, drop=FALSE]),
            InvS, TRUE)
        Fji[,i] <- -0.5*MhalDist - (p/2)*log(2*pi) - 0.5*logdetD
    }
    
    Fji <- sweep(Fji, 2, log(pivec), '+')
    Fjmax <- apply(Fji, 1, max)
    Fji <- sweep(Fji, 1, Fjmax, '-')
    Fji <- exp(Fji)
    tau <- sweep(Fji, 1, rowSums(Fji), '/')
    return(tau)
}

Mstep.mcfa <- function(Y, g, q, pivec, A, xi, omega, D, ...){
    p <- ncol(Y); if(is.null(p)) p <- 1
    n <- nrow(Y)
    tau <- tau.mcfa(Y, g, q, pivec, A, xi, omega, D)
    InvD <- diag(1/diag(D))
    A1 <- array(0, c(p,q))
    A2 <- array(0, c(q,q))
    Di <- array(0, c(p))
    for (i in seq_len(g)){
        gamma <- (InvD - InvD %*% A %*%
            chol.inv(chol.inv(omega[,,i]) + t(A) %*% InvD %*% A) %*%
            t(A) %*% InvD) %*% A %*% omega[,, i]
        ti <- sum(tau[,i])
        XI <- xi[, i, drop=FALSE]
        tY <- sweep(Y, 1, tau[,i], '*')
        Y_AXI <- sweep(Y, 2, A %*% XI, '-')
        tY_AXI <- sweep(Y_AXI, 1, tau[,i], '*')
        xi[, i] <- XI + t(gamma) %*% as.matrix(colSums(tY_AXI))/ti
        omega[,, i] <- (diag(q) - t(gamma) %*% A) %*% omega[,, i] +
            t(gamma) %*% t(Y_AXI) %*% tY_AXI %*% gamma/ti -
            (XI - xi[,i]) %*% t(XI - xi[,i])
        A1 <- A1 + colSums(tY) %*% t(XI) + t(Y) %*% tY_AXI %*% gamma
        A2 <- A2 + (omega[,, i] + xi[,i] %*% t(xi[,i]))* ti
        Di <- Di + diag(t(Y) %*% tY)
        pivec[i] <- ti/n
    }
    A <- try(A1 %*% chol.inv(A2))
    if(class(A)[1] == "try-error"){
        model <- "tried to invert a ill-conditioned or a singular matrix"
        class(model) <- "error"
        return(model)
    }
    D <- diag(Di - rowSums( (A %*% A2) * A ))/n
    model <- list(g=g, q=q, pivec=pivec, A=A, xi=xi, omega=omega, D=D)
    return(model)
}

init.est.para.mcfa <- function(Y, g, q, start, init.method="eigenA", ...){
    
    p <- ncol(Y)
    n <- nrow(Y)
    xi <- array(NA, c(q, g))
    omega <- array(NA, c(q, q, g))
    pivec <- array(NA, c(1, g))
    
    if(init.method=="randA")
    {
        A <- matrix(stats::rnorm(p*q), nrow=p, ncol=q)
        C <- chol(t(A) %*% A)
        A <- A %*% solve(C)
        D <- diag(diag(stats::cov(Y)))
        SqrtD <- diag(sqrt(diag(D)))
        InvSqrtD <- diag(1/diag(SqrtD))
        for(i in seq_len(g)) {
            indices  <- which(start == i)
            pivec[i] <- length(indices)/n
            uiT <- Y[indices,] %*% A
            xi[, i]  <- apply(uiT, 2, mean)
            Si <- stats::cov(Y[indices,])
            eig.list    <- try(eigen(InvSqrtD %*% Si %*% InvSqrtD), TRUE)
            H <- eig.list$vectors
            sort.lambda <- sort(eig.list$values, 
                                decreasing = TRUE, index.return=TRUE)
            lambda <- sort.lambda$x
            ix.lambda <- sort.lambda$ix
            if (q==p) {
                sigma2 <- 0
            }else {
                sigma2 <- mean(lambda[(q+1):p])
            }
            if(q==1) {
                omega[,, i] <- t(A) %*% SqrtD %*% H[, ix.lambda[seq_len(q)]] %*%
                    diag((lambda[seq_len(q)] - sigma2), q) %*%
                    t(H[, ix.lambda[seq_len(q)]]) %*% SqrtD %*% A
            }else{
                omega[,, i] <- t(A) %*% SqrtD %*% H[, ix.lambda[seq_len(q)]] %*%
                    diag((lambda[seq_len(q)] - sigma2))%*%
                    t(H[, ix.lambda[seq_len(q)]]) %*% SqrtD %*% A
            }
        }
    } else {
        Eig <- eigen(t(Y) %*% Y)
        A <- as.matrix(Eig$vectors[, seq_len(q)], ncol=q)
        for(i in seq_len(g)) {
            indices  <- which(start == i)
            pivec[i] <- length(indices)/n
            uiT <- Y[indices,] %*% A
            xi[, i]  <- apply(uiT, 2, mean)
            omega[,, i] <- stats::cov(t(t(A) %*% t(Y[indices,])) - xi[,i])
        }
        D <- diag(rep(mean(Eig$values[(q+1):p]), p))
    }
    model <- list(g=g, q=q, pivec=pivec, A=A, xi=xi, omega=omega, D=D)
    class(model) <- 'mcfa'
    return(model)
}

factor.scores.mcfa <- function(Y, g, q, pivec, A, xi, omega, D,
            tau=NULL, clust=NULL, ...){
    p <- ncol(Y); if(is.null(p)) p <- 1
    n <- nrow(Y)
    
    U <- array(0, c(n, q, g))
    gamma <- array(0, c(p, q, g))
    invD <- diag(1/diag(D))
    for (i in seq_len(g)){
        gamma[,, i] <- (invD - invD %*% A %*%
            chol.inv(chol.inv(omega[,, i]) + t(A) %*% invD %*% A) %*%
            t(A) %*% invD) %*% A %*% omega[,, i]
        U[,,i] <- t(replicate(n , xi[, i, drop=FALSE], 'matrix'))
        U[,,i] <- U[,,i] + 
            (Y - t(replicate(n, A%*% xi[, i, drop=FALSE], "matrix"))) %*%
            as.matrix(gamma[,, i])
    }
    
    if(is.null(tau))
        tau <- tau.mcfa(Y, g, q, pivec, A, xi, omega, D)
    if(is.null(clust))
        clust <- apply(tau, 1, which.max)
    UC <- array(0, c(n, q))
    Fmat <- array(0, c(n, q))
    for (i in seq_len(n)){
        UC[i, ] <- U[i,, clust[i]]
        Fmat[i, ] <- tau[i,] %*% t(matrix(U[i,,], c(q, g)))
    }
    return(list(U=U, UC=UC, Fmat=Fmat))
}


predict.mcfa <- function(object, x, ...){
    tau <- tau.mcfa(x, object$g, object$q, object$pivec, object$A, object$xi,
                    object$omega, object$D)
    clust <- apply(tau, 1, which.max)
    clust
}
