# function to compute efficient scores without family structure
geneStat <- function(x, Z=NULL, GenABEL.data=NULL, alpha, theta=NULL, gtime,
										 delta, beta=0, nCores=1,
										 FUN=function(Uij,weight){sum((colSums(Uij)*weight)^2)},
										 geneSet)
{
  if(sum(is.infinite(gtime)) >= 1)
	   ntps <- nlevels(as.factor(gtime)) - 1
	else
	   ntps <- nlevels(as.factor(gtime))
	
	# if user provided GenABEL object
    if (class(GenABEL.data) == "gwaa.data") {
      if (!is.null(x))
			  x <- as.matrix(as.numeric(GenABEL.data@gtdata[,x]))
		  else
			  x <- as.matrix(as.numeric(GenABEL.data@gtdata))
    }
  
  # transfer alphas into gammas and prevent numerical crash by setting zero to 10^-8
  alpha[which(alpha[1:ntps] < 10^-8)] <- 10^-8
  alpha[1:ntps] <- log(-log(alpha[1:ntps]))
  
	# make it as one vector for the nuisance parameter
	Params <- c(alpha, theta)

  # template results
  resStat <- list()
  resU <- list()
	# compute the gene level statistics for each SNPs set for genes
  if (is.null(colnames(x)))
    stop("SNP matrix does not have SNP ids, i.e., no colume names are provied")
  for (i in seq_len(length(geneSet))) {
    subSet_X <- as.matrix(x[, geneSet[[i]][[1]]], ncol=length(geneSet[[i]][[1]]) )
    # compute the efficient score by c++ function
    Uij <- .Call("_groupedSurv_effScore", PACKAGE = "groupedSurv", beta, Params, 
    subSet_X, Z, gtime, delta, ntps, nCores, TRUE)$uScore
    colnames(Uij) <- geneSet[[i]][[1]]
		resStat[[i]] <- FUN(Uij, geneSet[[i]][[2]]) 
    resU[[i]] <- Uij#, weight = geneSet[[i]][[2]])
  }
  
  res <- list()
	res$stat <- resStat
	res$U <- resU
  if(!is.null(names(geneSet)))
	{
		names(res$stat) <- names(geneSet)
		names(res$U) <- names(geneSet)
	}
	else
	{
		names(res$stat) <- paste0("gene ", 1:length(resStat))
		names(res$U) <- paste0("gene ", 1:length(resU))
	}
	return(res)
}







