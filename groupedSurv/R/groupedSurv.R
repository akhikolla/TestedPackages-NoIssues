 #function to compute efficient scores without family structure
.effScore <- function(x, Z=NULL, GenABEL.data=NULL, alpha, theta=NULL,
						      gtime, delta, beta=0, ntps=max(gtime), nCores=1, reScore=FALSE){

	# if user provided GenABEL object
  if (class(GenABEL.data) == "gwaa.data") {
    if(is.null(x))
			x <- as.matrix(as.numeric(GenABEL.data@gtdata))
		else
			x <- as.matrix(as.numeric(GenABEL.data@gtdata[,x]))
	}
  
	if(class(x) != "matrix"){
			x <- as.matrix(x, ncol=ncol(x))
	}
  
	if(class(Z) != "matrix"){
			Z <- as.matrix(Z,ncol=ncol(Z))
	}

  # transfer alphas into gammas
  alpha[which(alpha[1:ntps] < 10^-8)] <- 10^-8
  alpha[1:ntps] <- log(-log(alpha[1:ntps]))
  
	# make it as one vector for the nuisance parameter
	Params <- c(alpha, theta)

  # compute the efficient score by c++ function
  stat <- .Call("_groupedSurv_effScore", PACKAGE = "groupedSurv", beta, Params, 
    x, Z, gtime, delta, ntps, nCores, reScore)
 
  # compute p-values	
  pvalue <- pchisq(q = stat$sStatics, df = 1, lower.tail = FALSE)
  if(length(pvalue)>10){
      FDR <- qvalue(pvalue)$qvalues
      FWER <- p.adjust(pvalue, "bonferroni")
      res <- data.frame(stat = stat$sStatics, pvalue = pvalue, FDR=FDR, FWER=FWER)
	}else{
      cat("pi0 is can not be evaluated due to small number of p-values\nFDR and FWER is not calculated.\n")
      res <- data.frame(stat = stat$sStatics, pvalue = pvalue)}


	#res <- data.frame(stat = stat$sStatics, pvalue = pvalue, FDR=FDR, FWER=FWER)

  	
	# add column names
	rownames(res) <- colnames(x)
 
  # return results	
	if (!reScore) {
    return(res)
  } else {
    res <- list(stat = res, uScore = stat$uScore)
    return(res)
  }
}



groupedSurv <- function(x, Z=NULL, GenABEL.data=NULL, alpha, theta=NULL,
						      gtime, delta, beta=0, nCores=1, reScore=FALSE){
    if(sum(is.infinite(gtime)) >= 1)
      ntps <- nlevels(as.factor(gtime)) - 1
    else
      ntps <- nlevels(as.factor(gtime))
    
    if(class(Z) == "character")
      Z <- GenABEL.data@phdata[,Z]
   
	  .effScore(x, Z, GenABEL.data, alpha, theta,
							      gtime, delta, beta, ntps, nCores, reScore)
}



