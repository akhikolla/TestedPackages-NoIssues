## function to compute efficient score for data with family structure

.effScoreFam <- function(x, fam_group, fam_role, alpha, var, gtime, delta, beta=0, nCores=1){
	m <- length(unique(fam_group))

	registerDoParallel(cores=nCores)
  if(class(x) != "matrix"){
		SNPnumber <- 1
	  x <- as.matrix(x, ncol = 1)
	}else{
		SNPnumber <- ncol(x)
	}
  i <- 1  
	U <- foreach( i = 1:SNPnumber, .combine=rbind ) %dopar% {
	  res <- .Call("_groupedSurv_effScoreFam", PACKAGE = "groupedSurv", beta, var, 
      fam_group, alpha, gtime, delta, x[,i], fam_role, m)
    # print(res)
    U_beta <- res$beta_score
		#return(U_beta)
    C_Beta_Eta <- res$C_Beta_Eta
    V_Eta_Eta <- res$V_Eta_Eta
    U_Eta <- res$eta_score
    # print(U_Eta)
  	U_eff <- U_beta - C_Beta_Eta %*% solve(V_Eta_Eta) %*% U_Eta
  	#U_eff
	}
	t(U)
}

groupedSurvFam <- function(x, fam_group, fam_role, alpha, var, gtime, delta, beta=0, nCores=1){
	m <- length(unique(fam_group))
	.effScoreFam(x, fam_group, fam_role, alpha, var, gtime, delta, beta, nCores)
}


