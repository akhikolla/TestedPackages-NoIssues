predictions <- function(classif, x=matrix(0,nrow=1,ncol=1)){
	seed = get_seed()
	if(length(classif@icl)==0){
		warning('bosclassif returned an empty object, due to empty clusters. Predictions cannot be run.')
	}
	else{
		res <- prediction(classif, x, seed)
		return(res)
	}	
}