bosclassif <- function(x=matrix(0,nrow=1,ncol=1), y, idx_list=c(1), kr, kc=0, init, 
						nbSEM, nbSEMburn, nbindmini, m=0, percentRandomB=0){
	idx_list <- idx_list - 1 # patch for indexes
	checkParamsClassif(x, init, nbSEM, nbSEMburn)
	seed = get_seed()
	if(kc[1]!=0){
		res <- classif(xMat=x, y, myList=idx_list, kr, kc, init, nbSEM, nbSEMburn, 
			nbindmini, m=m, percentRandomB=percentRandomB, seed=seed)
	}
	else{
		res <- classifM(xMat=x, y, myList=idx_list, kr, init, nbSEM, nbSEMburn, 
			nbindmini, m=m, seed=seed)
	}
	if(length(res@icl)==0){
		warning('The algorithm found a spurious solution with empty clusters. You can: 1) Run the algorithm with another type of initialisation, 2) If you run the algorithm with init to "random" or "randomBurnin", running it again will change the initialisation , 3) Run the algorithm with a smaller argument kc.')
	}
	return(res)
}