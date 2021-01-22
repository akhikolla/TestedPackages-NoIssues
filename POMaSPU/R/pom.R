#' @useDynLib POMaSPU
#' @export POMaSPU
#' @importFrom Rcpp evalCpp
#' @importFrom MASS polr
#' @importFrom matrixStats rowMins
#' @importFrom stats fitted as.formula model.matrix






POMaSPU <- function(Y, Y.level, pheno, Z = NULL, pow = c(1:8, Inf), n.perm)
{
	if (n.perm < 1e7) a <- smallB(Y = Y, Y.level = Y.level, X = pheno, Z = Z,  pow = pow, B = n.perm) else
	a <- largeB(Y = Y, Y.level = Y.level, X = pheno, Z = Z,  pow = pow, B = n.perm)
	return(a)
}




rcum <- function(z){
  n.z <- length(z)
  rsum <- numeric(n.z - 1)
  for (i in 1:(n.z-1)) rsum[i] <- sum(z[i:(i+1)])
  rsum
}






smallB <- function(Y, Y.level, X, Z,  pow = c(1:8, Inf), B)
{
	# Null statistics generated using permutation
	
	n <- length(Y)
	Y <- factor(Y, levels = Y.level, ordered = TRUE)
	d <- table(Y)
	n.group <- length(d)
	indic <- data.frame(model.matrix(~Y-1)) # group indicator matrix
   
	if (is.null(Z)){
		r <- c(0, cumsum(as.vector(d)))/sum(d)
		U <- (1 - rcum(r)) %*% (t(indic) %*% X)
	} else {
		Z <- as.matrix(Z)
		colnames(Z) <- xnam <- paste0("z", 1:ncol(Z))
		fmla <- as.formula(paste("Y ~ ", paste(xnam, collapse= "+")))
		fit.logit <- polr(fmla, data = cbind.data.frame(Y, Z))
		fitted.prob <- fitted(fit.logit)
		r <- apply(fitted.prob, 1, function(z) c(0, cumsum(z))) 
		rsum0 <- apply(r, 2, function(w) 1- rcum(w))
		rsum <- apply(t(rsum0) * indic, 1, sum) 
		U <- crossprod(rsum, as.matrix(X))
	}
  
  
	# Observed statistics
	Ts <- rep(NA,length(pow))
	for (j in 1:length(pow)){
		if (pow[j] < Inf) Ts[j] <- sum(U^pow[j]) else Ts[j] <- max(abs(U))
	}
	
	# Permutation
	T0s <- matrix(0, nrow = B, ncol = length(pow))
	for (b in 1:B){
			  			  
 	        if (is.null(Z)){ 
		Y0 <- sample(Y, n)
		indic0 <- data.frame(model.matrix(~Y0 - 1)) 
		U0 <- (1 - rcum(r)) %*% (crossprod(indic0, X))
		} else {
		U0 <- crossprod(rsum, as.matrix(X[sample(1:n, n), ]))
		}
					
		# test stat's:
		for(j in 1:length(pow)){
			if (pow[j] < Inf) T0s[b, j] <- round(sum(U0^pow[j]), digits = 4)
			if (pow[j] == Inf) T0s[b, j] <- round(max(abs(U0)), digits = 4)
		}
	}
			
	pPerm0 <- rep(NA, length(pow));
	for(j in 1:length(pow)) {
		pPerm0[j] <- sum(abs(T0s[, j]) >= abs(Ts[j]))/B
	}
		
	P0s <- apply(T0s, 2, function(z) (B-rank(abs(z))+1)/B) 
	minP0s <- rowMins(P0s)
	Paspu <- (sum(minP0s <= min(pPerm0)) + 1)/(B + 1)
	
	namePerm0 <- paste("SPU", pow,  sep=".")
	pvs <- c(pPerm0, Paspu)
	names(pvs) <- c(namePerm0, "aSPU")
	return(pvs)
 }




largeB <- function(Y, Y.level, X, Z,  pow = c(1:8, Inf), B)
{
	pow_ <- ifelse(pow == Inf, 0, pow)
	n <- length(Y)
	X <- as.matrix(X) 
	Y <- factor(Y, levels = Y.level, ordered = TRUE)
	d <- table(Y)
	n.group <- length(d)
	indic <- data.frame(model.matrix(~Y-1)) # group indicator matrix
  
	Z <- as.matrix(Z)
	colnames(Z) <- xnam <- paste0("z", 1:ncol(Z))
	fmla <- as.formula(paste("Y ~ ", paste(xnam, collapse= "+")))
	fit.logit <- polr(fmla, data = cbind.data.frame(Y, Z))
	fitted.prob <- fitted(fit.logit)
	r <- apply(fitted.prob, 1, function(z) c(0, cumsum(z))) 
	rsum0 <- apply(r, 2, function(w) 1- rcum(w))
	rsum <- apply(t(rsum0) * indic, 1, sum) 
	U <- crossprod(rsum, X)
  
  
	##observed statistics
	Ts <-rep(NA,length(pow))
	for (j in 1:length(pow)){
		if (pow[j]<Inf) Ts[j]  <- sum(U^pow[j]) else Ts[j]  <- max(abs(U))
	}
	  
	P <- permhigh(as.matrix(Ts), as.matrix(rsum), as.matrix(X), as.matrix(pow_), B)
	pvs <- c(P$P_SPU, P$P_aSPU)
	namePerm0 <- paste("SPU", pow,  sep=".")
	names(pvs) <- c(namePerm0, "aSPU")
	return(pvs)
}

