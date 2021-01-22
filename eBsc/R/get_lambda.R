get_lambda <- function(zero_range, coefs, eigens, rhos, neBsc, q){
    eps <- .Machine$double.eps
    aux <- try(uniroot(Elambda, zero_range, coefs = coefs, eigens = eigens, rhos = rhos, neBsc = neBsc, q = q, tol = eps, maxiter = 10^5, extendInt = "yes")$root)
    if(class(aux)=="try-error"){
        out <- 5; message("Warning: non-optimal smoothing parameter.")
    }else{
        out <- aux}
    out
}

