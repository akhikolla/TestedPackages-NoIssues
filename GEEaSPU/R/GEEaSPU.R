#' @useDynLib GEEaSPU
#' @export GEEaSPU GEEaSPUset GEEaSPUpath
#' @importFrom Rcpp evalCpp
#' @importFrom gee gee
#' @importFrom stats cov fitted formula


gin <- function (X, tol = sqrt(.Machine$double.eps)) {

    if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
        stop("'X' must be a numeric or complex matrix")
    if (!is.matrix(X)) 
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X)) 
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    if (all(Positive)) 
        Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive)) 
        array(0, dim(X)[2L:1L])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
        t(Xsvd$u[, Positive, drop = FALSE]))
}




longForm <-function(dat, k){
  
  p <- ifelse(is.vector(dat), 1, ncol(dat))
  n <- ifelse(is.vector(dat), length(dat), nrow(dat))
  newdat <- matrix(0, ncol = p*k, nrow = n*k)
  for (r in 1:k) newdat[seq(r, n*k, k),(p*(r-1)+1):(p*r)] <- as.matrix(dat)
  return(newdat)
}






sB <- function(pheno, geno, Z , model, corstr, pow, pow2, n.perm, null.type, score.test){
	  
	if (model=="binomial" && null.type=="perm") stop("Please use a simulation method for binary phenotypes.")     
	  
	p <- ncol(geno); geno <- as.matrix(geno)
	pow_ <- as.matrix(ifelse(pow==Inf, 0, pow))
	pow2_ <- as.matrix(ifelse(pow2==Inf, 0, pow2))

	if (is.null(dim(pheno))){
		n <- length(pheno); k <- 1 
	} else {
		n <- nrow(pheno); k <- ncol(pheno)
	}
	id <- rep(1:n, each = k) 
	longX <- longForm(geno, k);
	longy <- matrix(c( t(pheno)), n*k, 1)    
	if (is.null(Z)) {
		longZ <- longForm(matrix(1, ncol = 1, nrow = n), k); mf <- formula(longy ~ 0 + longZ)
	} else { 
		longZ <- longForm(Z, k); mf <- formula(longy ~ longZ)
	}
	x <- cbind(longZ, longX)

	geefit <- gee(mf, id = id, family = model, corstr = corstr)
	R <- geefit$working.correlation
	muy <- fitted(geefit)
	res <- matrix(geefit$residuals, ncol = k, byrow = T)
	covres <- cov(res,use = "pairwise.complete.obs")
	if (model == "gaussian"){	
		v0 <- rep(1,n*k); f=0
	} else {
		v0 <- muy*(1-muy); f=1 
	}

	if (null.type=="perm" && score.test==FALSE)
		fit <- perm(gin(R), geno, res, n, k, p, pow_, pow2_, n.perm)	    
	if (null.type=="perm" && score.test)
		fit <- perm_score(gin(R), geno, x, res, covres, n,k,p,ncol(x), pow_, pow2_, n.perm) 
	if (null.type=="sim" && score.test==FALSE)
		fit <-  sim(f, v0, gin(R), geno, x, res, covres, n,k,p,ncol(x), pow_, pow2_, n.perm)
	if (null.type=="sim" && score.test)
		fit <- sim_score(f, v0, gin(R), geno, x, res, covres, n,k,p,ncol(x), pow_, pow2_, n.perm)

	return(fit)
}
  
  
hB <- function(pheno, geno, Z , model, corstr, pow, pow2, n.perm, null.type, score.test){
	  
	if (model=="binomial" && null.type=="perm")  stop("Please use a simulation method for binary phenotypes.")     
	  
	p <- ncol(geno); geno <- as.matrix(geno)
	pow_ <- as.matrix(ifelse(pow == Inf, 0, pow))
	pow2_ <- as.matrix(ifelse(pow2 == Inf, 0, pow2))

	if (is.null(dim(pheno))){
		n <- length(pheno); k <- 1 
	} else {
		n <- nrow(pheno); k <- ncol(pheno)
	}
	id <- rep(1:n, each = k) 
	longX <- longForm(geno, k);
	longy <- matrix(c( t(pheno)), n*k, 1)    
	if (is.null(Z)) {
		longZ <- longForm(matrix(1,ncol=1, nrow=n), k); mf <- formula(longy ~ 0 + longZ)
	  } else { 
		longZ <- longForm(Z, k); mf <- formula(longy ~ longZ)
	}
	x <- cbind(longZ, longX)

	geefit <- gee(mf, id = id, family = model, corstr = corstr)
	R <- geefit$working.correlation
	muy <- fitted(geefit)
	res <- matrix(geefit$residuals, ncol = k, byrow = T)
	covres <- cov(res,use = "pairwise.complete.obs")
	if (model =="gaussian"){	
		v0 <- rep(1,n*k); f = 0
	} else {
		v0 <- muy*(1-muy); f = 1 
	}

	if (null.type == "perm" && score.test == FALSE)
		fit <- permhigh(gin(R), geno, res, n, k, p, pow_, pow2_, n.perm)	    
	if (null.type=="perm" && score.test)
		fit <- permhigh_score(gin(R), geno, x, res, covres, n,k,p,ncol(x), pow_, pow2_, n.perm) 
	if (null.type=="sim" && score.test==FALSE)
		fit <-  simhigh(f, v0, gin(R), geno, x, res, covres, n,k,p,ncol(x), pow_, pow2_, n.perm)
	if (null.type=="sim" && score.test)
		fit <- simhigh_score(f, v0, gin(R), geno, x, res, covres, n,k,p,ncol(x), pow_, pow2_, n.perm)

	return(fit)
}
 


pathB <- function(pheno, geno, nSNPs, Z , model, corstr, pow, pow2, pow3, n.perm){

	nSNPs <- nSNPs[nSNPs>=1]
	nGenes <- length(nSNPs)
	p <- ncol(geno); geno <- as.matrix(geno)
	pow_ <- as.matrix(ifelse(pow==Inf, 0, pow))
	pow2_ <- as.matrix(ifelse(pow2==Inf, 0, pow2))
	pow3_ <- as.matrix(ifelse(pow3==Inf, 0, pow3))

	if (is.null(dim(pheno))){
		n <- length(pheno); k <- 1 
	} else {
		n <- nrow(pheno); k <- ncol(pheno)
	}
	id <-rep(1:n, each = k) 
	longy <- matrix(c( t(pheno)), n*k, 1)    
	if (is.null(Z)) {
		longZ <- longForm(matrix(1,ncol=1, nrow=n), k); mf <- formula(longy ~ 0 + longZ)
	} else { 
		longZ <- longForm(Z, k); mf <- formula(longy ~ longZ)
	}
	geefit <- gee(mf, id = id, family = model, corstr = corstr)
	R <- geefit$working.correlation
	muy <- fitted(geefit)
	res <- matrix( geefit$residuals, ncol = k, byrow = T)
	fit <- permpath(gin(R), geno, res, nSNPs, nGenes, n, k, p, pow_, pow2_,  pow3_, n.perm)
return(fit)
}
 











GEEaSPU <-function(pheno, geno, Z = NULL, model = "gaussian", corstr = "independence", pow = c(1:8, Inf), n.perm = 1000, null.type = "perm", score.test = FALSE)
{
	geno <- matrix(geno, ncol=1)
	if(n.perm <=5e7){
		fit <- sB(pheno = pheno, geno = geno, Z = Z, model = model, corstr = corstr, pow=1, pow2=pow, n.perm = n.perm, null.type = null.type, score.test = score.test)
	  } else{
		fit <- hB(pheno = pheno, geno = geno, Z = Z, model = model, corstr = corstr, pow=1, pow2=pow, n.perm = n.perm, null.type = null.type, score.test = score.test)
	}
	if (score.test == TRUE){
		pvs <- c(fit$P_SPU, fit$P_aSPU, fit$P_aSPU.score)
		pvs[pvs==0] <- 1/n.perm
		names(pvs) <- c(paste("SPU", pow, sep = "."), "Score", "aSPU", "aSPU_score")
	} else {
		pvs <- c(fit$P_SPU, fit$P_aSPU)
		pvs[pvs==0] <- 1/n.perm
		names(pvs) <- c(paste("SPU", pow, sep = "."), "aSPU")
	}
	return(pvs)
}


GEEaSPUset <-function(pheno, geno, Z = NULL, model = "gaussian", corstr = "independence", pow = c(1,2,4,8), pow2 = c(1,2,4,8), n.perm = 1000, null.type = "perm", score.test = FALSE)
{	  
	if(n.perm <=1e7){
		fit <- sB(pheno = pheno, geno = geno, Z = Z, model = model, corstr = corstr, pow = pow, pow2 = pow2, n.perm = n.perm, null.type = null.type, score.test = score.test)
	} else{
		fit <- hB(pheno = pheno, geno = geno, Z = Z, model = model, corstr = corstr, pow = pow, pow2 = pow2, n.perm = n.perm, null.type = null.type, score.test = score.test)
	}
	if (score.test == TRUE){
		pvs <- c(fit$P_SPU, fit$P_aSPU, fit$P_aSPU.score)
		pvs[pvs==0] <- 1/n.perm
		names(pvs) <- c(sapply(pow, function(z) paste("SPU", z, pow2, sep = ".")), "Score", "aSPUset", "aSPUset_score")
	} else {
		pvs <- c(fit$P_SPU, fit$P_aSPU)
		pvs[pvs==0] <- 1/n.perm
		names(pvs) <- c(sapply(pow, function(z) paste("SPU", z, pow2, sep = ".")), "aSPUset")
	}
	return(pvs)

}



GEEaSPUpath <-function(pheno, geno, nSNPs, Z = NULL, corstr = "independence",  pow = c(1,2,4,8), pow2 = c(1,2,4,8), pow3 = c(1,2,4,8), n.perm = 1000)
{	
  model <- "gaussian" ; 	  
  fit <- pathB(pheno = pheno, geno = geno, nSNPs = nSNPs, Z = Z, model = model, corstr = corstr, pow = pow, pow2 = pow2, pow3 = pow3, n.perm = n.perm)
	pvs <- c(fit$P_SPU, fit$P_aSPUpath)
	pvs[pvs==0] <- 1/n.perm
	nam <- c(sapply(pow, function(x) paste(x, pow2, sep = ".")))
	names(pvs) <- c(sapply(nam, function(x) paste("SPU", x, pow3, sep = ".")), "aSPUpath")	
	return(pvs)
}


