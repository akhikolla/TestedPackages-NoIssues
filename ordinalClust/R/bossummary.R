
bossummary <- function(object){
	
	m <- object@m

	D <- length(m)


	nb.V <- length(unique(object@zr))
	gammas <- rep(0,nb.V)
	for(i in 1:nb.V){
		gammas[i] <- length(which(object@zr==i))/length(object@zr)
	}

	

	
	if(object@name == "ClassifM"){
		print("############################ CLASSIFICATION RESULTS ############################")
		
	}
	if(object@name == "Classif"){
		print("############################ CLASSIFICATION RESULTS ############################")
		print("########################## Columns mixing proportions ##########################")
		for(id in 1:D){
			print(paste("* group", id, "*"))
			nb.W <- length(unique(object@zc[[id]]))

			rhos <- c()
			if(nb.W!=1){
				for(i in 1:(nb.W)){
					rho <- length(which(object@zc[[id]]==i))/length(object@zc[[id]])
					rhos <- c(rhos,rho)
				}
			}
			print(rhos)
		}
	}
	if(object@name == "Coclust"){
		print("############################ CO-CLUSTERING RESULTS #############################")
		print("########################## Columns mixing proportions ##########################")
		for(id in 1:D){
			print(paste("* group", id, "*"))
			nb.W <- length(unique(object@zc[[id]]))

			rhos <- c()
			if(nb.W!=1){
				for(i in 1:(nb.W)){
					rho <- length(which(object@zc[[id]]==i))/length(object@zc[[id]])
					rhos <- c(rhos,rho)
				}
			}
			print(rhos)
		}
	}
	if(object@name == "Clust"){
		print("############################## CLUSTERING RESULTS ##############################")
	}

	print("############################ Rows mixing proportions ###########################")
	print(gammas)

	print("################################ BOS parameters ################################")
	for(id in 1:D){
		print(paste("################################### Group", id, "####################################"))
		print(object@params[[id]])
	}


	
}