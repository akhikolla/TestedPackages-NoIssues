boscoclust <- function(x=matrix(0,nrow=1,ncol=1), idx_list=c(1), 
	kr, kc, init, nbSEM, nbSEMburn, 
	nbRepeat=1, nbindmini, m=0, percentRandomB=0){
	idx_list <- idx_list - 1 # patch for indexes
	checkParamsCoclust(x, init, nbSEM, nbSEMburn)
	seed = get_seed()
	res <- coclust(xMat=x, myList=idx_list, kr, kc, init, nbSEM, nbSEMburn, 
		nbRepeat, nbindmini, m=m, percentRandomB=percentRandomB, seed=seed)
	if(length(res@icl)==0){
		warning('The algorithm found a spurious solution with empty clusters. You can: 1) Run the algorithm with another type of initialisation, 2) If you run the algorithm with init to "random" or "randomBurnin", running it again will change the initialisation , 3) Run the algorithm with a smaller argument kr and/or kc.')
	}
	return(res)
}