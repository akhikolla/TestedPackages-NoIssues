bosclust <- function(x=matrix(0,nrow=1,ncol=1), idx_list=c(1), kr, init, 
					nbSEM, nbSEMburn, nbindmini, m=0, percentRandomB=0){
	idx_list <- idx_list - 1 # patch for indexes
	checkParamsClust(x, init, nbSEM, nbSEMburn)
	seed = get_seed()
	res <- clust(xMat=x, myList=idx_list, kr, init, nbSEM, nbSEMburn, nbindmini, 
				m=m, percentRandomB=percentRandomB, seed=seed)
	if(length(res@icl)==0){
		warning('The algorithm found a spurious solution with empty clusters. You can: 1) Run the algorithm with another type of initialisation, 2) If you run the algorithm with init to "random" or "randomBurnin", running it again will change the initialisation , 3) Run the algorithm with a smaller argument kr.')
	}
	return(res)
}