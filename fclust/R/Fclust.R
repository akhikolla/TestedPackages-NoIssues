Fclust <- function (X, k, type, ent, noise, stand, distance)
	{
	if (missing(X))
	stop("The data set must be given")
	if (is.null(X))
	stop("The data set X is empty")
	X=as.matrix(X)
	if (any(is.na(X)))
		stop("The data set X must not contain NA values")
	if (!is.numeric(X)) 
		stop("The data set X is not a numeric data.frame or matrix")
	if (missing(distance))
			distance=FALSE
		if (any(is.na(match(distance,c(FALSE,TRUE)))))
			{
			distance=FALSE
			cat("The default value distance=FALSE has been set ",fill=TRUE)
		}
	if (distance==FALSE)
		{	
		n=nrow(X)
		p=ncol(X)
		if (is.null(rownames(X)))
			rn=paste("Obj",1:n,sep=" ")
		else
			rn=rownames(X)
		if (is.null(colnames(X)))
			cn=paste("Var",1:p,sep=" ")
		else
			cn=colnames(X)
		if (missing(type))
			type="standard"
		if (any(is.na(match(type,c("standard","polynomial","gk","gkb","medoids")))))
			{
			type="standard"
			cat("The default value type=\"standard\" has been set ",fill=TRUE)
		}
		if (missing(ent))
			ent=FALSE
		if (any(is.na(match(ent,c(FALSE,TRUE)))))
			{
			ent=FALSE
			cat("The default value ent=FALSE has been set ",fill=TRUE)
		}
		if (missing(noise))
			noise=FALSE
		if (any(is.na(match(noise,c(FALSE,TRUE)))))
			{
			noise=FALSE
			cat("The default value noise=FALSE has been set ",fill=TRUE)
		}
		if (missing(k))
			{
			k=2
			cat("The number of clusters k is missing: k=2 will be used ",fill=TRUE)
		}
		if (!is.numeric(k)) 
			{
			k=2
			cat("The number of clusters k is not numeric: k=2 will be used ",fill=TRUE)
		}
		if ((k>ceiling(n/2)) || (k<2))
			{
			k=2
			cat("The number of clusters k must be an integer in {2, 3, ..., ceiling(n/2)}: k=2 will be used ",fill=TRUE)
		}
		if (k%%ceiling(k)>0)  
			{
			k=2
			cat("The number of clusters k must be an integer in {2, 3, ..., ceiling(n/2)}: k=2 will be used ",fill=TRUE)
		}
		if (missing(stand))
			stand=0
		if (!is.numeric(stand))
			stand=0
		if (stand==1)
			X=scale(X,center=TRUE,scale=TRUE)[,]

		if (noise==FALSE)
			{
			if ((type=="standard") & (ent==FALSE))
				{
				clust=FKM(X,k)
				cat("The standard FkM algorithm has been chosen", fill=TRUE)
				cat("The default options have been set, to specify different options, use FKM ", fill=TRUE)
			}
			if ((type=="standard") & (ent==TRUE))
				{
				clust=FKM.ent(X,k)
				cat("The FkM algorithm with entropy regularization has been chosen", fill=TRUE)
				cat("The default options have been set, to specify different options, use FKM.ent ", fill=TRUE)
			}

		  
		  
		  
		  
		  if ((type=="gk") & (ent==FALSE))
				{
				clust=FKM.gk(X,k)
				cat("The Gustafson and Kessel extension of the FkM algorithm has been chosen", fill=TRUE)
				cat("The default options have been set, to specify different options, use FKM.gk ", fill=TRUE)
			}
			if ((type=="gk") & (ent==TRUE))
				{
				clust=FKM.gk.ent(X,k)
				cat("The Gustafson and Kessel extension of the FkM algorithm with entropy regularization has been chosen", fill=TRUE)
				cat("The default options have been set, to specify different options, use FKM.gk.ent ", fill=TRUE)
			}
			if ((type=="gkb") & (ent==FALSE))
				{
				clust=FKM.gkb(X,k)
				cat("The Gustafson, Kessel and Babuska extension of the FkM algorithm has been chosen", fill=TRUE)
				cat("The default options have been set, to specify different options, use FKM.gkb ", fill=TRUE)
			}
			if ((type=="gkb") & (ent==TRUE))
				{
				clust=FKM.gkb.ent(X,k)
				cat("The Gustafson, Kessel and Babuska extension of the FkM algorithm with entropy regularization has been chosen", fill=TRUE)
				cat("The default options have been set, to specify different options, use FKM.gkb.ent ", fill=TRUE)
			}
			if (type=="polynomial") 
				{	
				if (ent==TRUE) cat("When type=\"polynomial\", ent=FALSE has been set ",fill=TRUE)
				clust=FKM.pf(X,k)
				cat("The FkM algorithm with polynomial fuzzifier has been chosen", fill=TRUE)
				cat("The default options have been set, to specify different options, use FKM.pf ", fill=TRUE)
			}
			if (type=="medoids") 
				{
				if (ent==TRUE) cat("When type=\"medoids\", ent=FALSE has been set ",fill=TRUE)
				clust=FKM.med(X,k)
				cat("The fuzzy k-medoids algorithm has been chosen", fill=TRUE)
				cat("The default options have been set, to specify different options, use FKM.med ", fill=TRUE)
			}
		}
		else{
			if ((type=="standard") & (ent==FALSE))
				{
				clust=FKM.noise(X,k)
				cat("The standard FkM algorithm with noise cluster has been chosen", fill=TRUE)
				cat("The default options have been set, to specify different options, use FKM.noise ", fill=TRUE)
			}
			if ((type=="standard") & (ent==TRUE))
				{
				clust=FKM.ent.noise(X,k)
				cat("The FkM algorithm with entropy regularization and noise cluster has been chosen", fill=TRUE)
				cat("The default options have been set, to specify different options, use FKM.ent.noise ", fill=TRUE)
			}
			if ((type=="gk") & (ent==FALSE))
				{
				clust=FKM.gk.noise(X,k)
				cat("The Gustafson and Kessel extension of the FkM algorithm with with noise cluster has been chosen", fill=TRUE)
				cat("The default options have been set, to specify different options, use FKM.gk.noise ", fill=TRUE)
			}
			if ((type=="gk") & (ent==TRUE))
				{
				clust=FKM.gk.ent.noise(X,k)
				cat("The Gustafson and Kessel extension of the FkM algorithm with entropy regularization and noise cluster has been chosen", fill=TRUE)
				cat("The default options have been set, to specify different options, use FKM.gk.ent.noise ", fill=TRUE)
			}
			if ((type=="gkb") & (ent==FALSE))
				{
				clust=FKM.gkb.noise(X,k)
				cat("The Gustafson, Kessel and Babuska extension of the FkM algorithm with noise cluster has been chosen", fill=TRUE)
				cat("The default options have been set, to specify different options, use FKM.gkb.noise ", fill=TRUE)
			}
			if ((type=="gkb") & (ent==TRUE))
				{
				clust=FKM.gkb.ent.noise(X,k)
				cat("The Gustafson, Kessel and Babuska extension of the FkM algorithm with entropy regularization and noise cluster has been chosen", fill=TRUE)
				cat("The default options have been set, to specify different options, use FKM.gkb.ent.noise ", fill=TRUE)
			}
			if (type=="polynomial")
				{
				if (ent==TRUE) 
					cat("When type=\"polynomial\", ent=FALSE has been set ", fill=TRUE)
				clust=FKM.pf.noise(X,k)
				cat("The FkM algorithm with polynomial fuzzifier and noise cluster has been chosen", fill=TRUE)
				cat("The default options have been set, to specify different options, use FKM.pf.noise ", fill=TRUE)
			}
			if (type=="medoids") 
				{
				if (ent==TRUE) 
					cat("When type=\"medoids\", ent=FALSE has been set ", fill=TRUE)
				clust=FKM.med.noise(X,k)
				cat("The fuzzy k-medoids algorithm with noise cluster has been chosen", fill=TRUE)
				cat("The default options have been set, to specify different options, use FKM.med.noise ", fill=TRUE)
			}
		}
	}
	else{
	  n=nrow(X)
	  if (n!=ncol(X))
	    stop("The data set X is not a a distance/dissimilarity matrix matrix")
	  if (sum(t(X)-X)>1e-6) 
	    stop("The data set X is not a a distance/dissimilarity matrix matrix")
	  if (sum(diag(X))>1e-6) 
	    stop("The data set X is not a a distance/dissimilarity matrix matrix")
	  if (is.null(rownames(X)))
			rn=paste("Obj",1:n,sep=" ")
		else
			rn=rownames(X)
		cn=rn
		
		if (missing(type))
		  type="standard"
		if (type!="standard")
		  cat("When distance=TRUE, type=\"standard\" has been set ", fill=TRUE)
		if (missing(ent))
		  ent=FALSE
		if (ent!=FALSE)
		  cat("When distance=TRUE, ent=FALSE has been set ", fill=TRUE)
		if (missing(stand))
		  stand=0
		if (stand!=0)
		  cat("When distance=TRUE, stand=0 has been set ", fill=TRUE)
		if (missing(noise))
		  noise=FALSE
		if (any(is.na(match(noise,c(FALSE,TRUE)))))
		{
		  noise=FALSE
		  cat("The default value noise=FALSE has been set ",fill=TRUE)
		}
		if (missing(k))
		{
		  k=2
		  cat("The number of clusters k is missing: k=2 will be used ",fill=TRUE)
		}
		if (!is.numeric(k)) 
		{
		  k=2
		  cat("The number of clusters k is not numeric: k=2 will be used ",fill=TRUE)
		}
		if ((k>ceiling(n/2)) || (k<2))
		{
		  k=2
		  cat("The number of clusters k must be an integer in {2, 3, ..., ceiling(n/2)}: k=2 will be used ",fill=TRUE)
		}
		if (k%%ceiling(k)>0)  
		{
		  k=2
		  cat("The number of clusters k must be an integer in {2, 3, ..., ceiling(n/2)}: k=2 will be used ",fill=TRUE)
		}
		if (noise==FALSE)
			{
			clust=NEFRC(D=X,k)
			cat("The non-euclidean fuzzy relational algorithm has been chosen", fill=TRUE)
			cat("The default options have been set, to specify different options, use NEFRC ", fill=TRUE)
		}
		else{
			clust=NEFRC.noise(D=X,k)
			cat("The non-euclidean fuzzy relational algorithm with noise cluster has been chosen", fill=TRUE)
			cat("The default options have been set, to specify different options, use NEFRC.noise ", fill=TRUE)
		}
	}
	return(clust)
}