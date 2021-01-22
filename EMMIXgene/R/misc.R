ari <- function (cls, hat.cls)
{
    if(length(cls) != length(hat.cls))
        stop('length of the two arguments should be equal')
    tab <- table(cls, hat.cls)
    if (sum(diag(tab)) == length(cls))
        return(1)
    A <- sum(choose(tab, 2))
    B <- sum(choose(rowSums(tab), 2))
    C <- sum(choose(colSums(tab), 2))
    D <- choose(sum(tab), 2)
    
    ARI <- (A - B*C/D) / ( (B+C)/2 - B*C/D)
    
    return(ARI)
}

# permutations
perm <- function (n, r, v = seq_len(n)){
    if (r == 1)
        matrix(v, n, 1)
    else if (n == 1)
        matrix(v, 1, r)
    else {
        X <- NULL
        for (i in seq_len(n))
            X <- rbind(X, cbind(v[i], perm(n - 1, r - 1, v[-i])))
        X
    }
}

# the minimum number of misallocations
err <- function(cls, hat.cls) {
    if(length(cls) != length(hat.cls))
        stop('length of the two arguments should be equal')
    tcls = rep(0, length(cls))
    labs = unique(cls)
    for(j in seq_len(length(labs)))
        tcls[cls==labs[j]] = j
    cls = tcls
    
    tcls = rep(0, length(hat.cls))
    labs = unique(hat.cls)
    for(j in seq_len(length(labs)))
        tcls[hat.cls==labs[j]] = j
    hat.cls = tcls
    
    labs = unique(c(hat.cls, cls))
    g = length(labs)
    x = perm(g, g)
    min.err = Inf
    for( j in seq_len(nrow(x) )) {
        newCls = rep(0, length(cls))
        for(k in seq_len(g))
            newCls[cls==labs[k]] = x[j, k]
        e = sum(newCls != hat.cls)
        if(e < min.err)
            min.err = e
    }
    return(min.err)
}

Starts <- function(Y, g, initClust, nkmeans, nrandom){
    n <- nrow(Y)
    starts <- NULL
    
    if(!is.null(initClust))
    {
        if( (is.null(dim(initClust))) && (length(initClust) != n) )
            stop('initClust has not specified correctly.')
        
        if( (!is.null(dim(initClust))) && (dim(initClust)[1] != n) )
            stop('initClust has not specified correctly.')
    }
    
    if(is.numeric(nkmeans) & ! is.null(nkmeans) & (nkmeans > 0))
    {
        for (i in seq_len(nkmeans))
            starts <- cbind(starts, stats::kmeans(Y, g)$cluster)
    }
    
    if(is.numeric(nrandom) & ! is.null(nrandom) & (nrandom > 0))
        starts <- cbind(starts, floor(g * matrix(stats::runif(n*nrandom), 
            c(n, nrandom))) + 1)
    
    starts <- cbind(initClust, starts)
    starts <- as.matrix(starts)
    
    if(is.null(starts))
        stop("Alteast one of Y, g, initClust, nkmeans, 
        nrandom has not been specified correctly.")
    
    return(starts)
}

chol.inv <- function(x, ...){
    C <- chol(x)
    inv.x <- chol2inv(C)
    return(inv.x)
}
