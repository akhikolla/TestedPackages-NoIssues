x.cl <- function(m, arl0, alpha=0.10, beta=0.05, H=200, A=1.5, B=50, 
                 Ninit=1000, Nfinal=30000) {
    if (arl0<=2) stop("arl0 must be greater than 2")
    Linf <- qnorm(1-1/(2*arl0)) 
    .cautious.design(list(chart="X"), m, arl0, alpha, beta, 
                     Linf, A, B, H, Ninit, Nfinal)
}

ewma.cl <- function(lambda, m, arl0, alpha=0.10, beta=0.05, H=200, A=1.5, B=50, 
                    Ninit=1000, Nfinal=30000) {
    lambda <- as.numeric(lambda)
    if ( (length(lambda) != 1) || (lambda <= 0) || (lambda >= 1) ) 
        stop("lambda must be in (0,1)")
    if (arl0<=2) stop("arl0 must be greater than 2")
    Linf <- spc::xewma.crit(lambda, arl0, sided="two") 
    .cautious.design(list(chart="EWMA", lambda=lambda), 
                     m, arl0, alpha, beta, 
                     Linf, A, B, H, Ninit, Nfinal)
}

cusum.cl <- function(k, m, arl0, alpha=0.10, beta=0.05, H=200, A=1.5, B=50, 
                     Ninit=1000, Nfinal=30000) {
    k <- as.numeric(k)
    if ( (length(k) != 1) || (k <= 0) ) 
        stop("k must be greater than 0")
    if (arl0<=2) stop("arl0 must be greater than 2")
    Linf <- spc::xcusum.crit(k, arl0, sided="two") 
    .cautious.design(list(chart="CUSUM", k=k), 
                     m, arl0, alpha, beta, 
                     Linf, A, B, H, Ninit, Nfinal)
}


.cautious.design <- function(chart, m, arl0, alpha, beta, 
                             Linf, A, B, H, Ninit, Nfinal) {
    H <- round(H)
    Ninit <- round(Ninit)
    Nfinal <- round(Nfinal)
    if (m<10) stop("at least 10 Phase I observations must be available") 
    if ((alpha<=0) || (alpha>=1)) stop("alpha must be in (0,1)")  
    if ((beta<=0) || (beta>=1)) stop("beta must be in (0,1)")
    if (H<10) stop("H must be al least 10")
    if (Ninit<10) stop("Ninit must be at least 10")
    if (Nfinal<100) stop("Nfinal must be at least 100")
    mkChart(chart, m,  A, B, arl0, Linf, alpha, beta, H, Ninit, Nfinal)
}

