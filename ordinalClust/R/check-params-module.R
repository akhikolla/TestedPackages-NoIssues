checkParamsClust <- function(x, init, nbSEM, nbSEMburn){
	##### generic things #######
	generical(x, init, nbSEM, nbSEMburn)
}

checkParamsCoclust <- function(x, init, nbSEM, nbSEMburn){
	##### generic things #######
	generical(x, init, nbSEM, nbSEMburn)
}

checkParamsClassif <- function(x, init, nbSEM, nbSEMburn){
	##### generic things #######
	generical(x, init, nbSEM, nbSEMburn)
	#if((nrow(x)!=length(y))||(dim(functionalData)[1]!=length(y))) stop('y must have nrow(x) as number of elements.')
}




######################## UTILS ######################
generical <- function(x, init, nbSEM, nbSEMburn){
	if(!any(init=="kmeans") & !any(init=="random") & !any(init=="randomBurnin")) stop('init argument should be "kmeans", "random" or "randomBurnin".')
	if(nbSEM<=nbSEMburn) stop('nbSEMburn must be inferior to nbSEM.')
}

