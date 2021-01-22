bosplot1 <- function(object){

	m <- object@m

	nb.V <- length(unique(object@zr))
	gammas <- rep(0,nb.V)
	for(i in 1:nb.V){
		gammas[i] <- length(which(object@zr==i))/length(object@zr)
	}


	
	if(object@name == "ClassifM"){
		D <- length(m)
		par(mfrow=c(1,D))
		for(id in 1:D){
			tmp <- object@xhat[[id]][sort(object@zr,index.return=TRUE)$ix,1:ncol(object@xhat[[id]])]
			par(xpd=TRUE, mar=c(2,1,4,7))
			image(t(tmp),xaxt='n',yaxt='n', main='classification',cex.main=1.5,col = gray(m[id]:1/m[id]))


			legend(x=1.01,y=0.65,legend = 1:m[id],col = gray(m[id]:1/m[id]),pch=15,bty="n",cex=1)
			par(xpd=FALSE)
			sum.gammas <- 0
			for(i in 1:(nb.V-1)){
				sum.gammas <- sum.gammas + gammas[i]
				abline(h=sum.gammas,lwd=3, col="red")
			}
		}
	}
	if(object@name == "Classif"){
		D <- length(m)
		par(mfrow=c(1,D))
		for(id in 1:D){
			tmp <- object@xhat[[id]][sort(object@zr,index.return=TRUE)$ix,
					sort(object@zc[[id]],index.return=TRUE)$ix]
			par(xpd=TRUE, mar=c(2,1,4,7))
			image(t(tmp),xaxt='n',yaxt='n', main='classification',cex.main=1.5,col = gray(m[id]:1/m[id]))


			legend(x=1.01,y=0.65,legend = 1:m[id],col = gray(m[id]:1/m[id]),pch=15,bty="n",cex=1)
			par(xpd=FALSE)
			nb.W <- length(unique(object@zc[[id]]))

			sum.rho <- 0
			if(nb.W!=1){
				for(i in 1:(nb.W-1)){
					rho <- length(which(object@zc[[id]]==i))/length(object@zc[[id]])
					sum.rho <- sum.rho + rho
					abline(v=sum.rho,lwd=3, col="red")
				}
			}
			

			sum.gammas <- 0
			for(i in 1:(nb.V-1)){
				sum.gammas <- sum.gammas + gammas[i]
				abline(h=sum.gammas,lwd=3, col="red")
			}
		}
	}
	if(object@name == "Coclust"){
		D <- length(m)
		par(mfrow=c(1,D))
		for(id in 1:D){
			tmp <- object@xhat[[id]][sort(object@zr,index.return=TRUE)$ix,
					sort(object@zc[[id]],index.return=TRUE)$ix]
			par(xpd=TRUE, mar=c(2,1,4,7))
			image(t(tmp),xaxt='n',yaxt='n', main='co-clustering',cex.main=1.5,col = gray(m[id]:1/m[id]))


			legend(x=1.01,y=0.65,legend = 1:m[id],col = gray(m[id]:1/m[id]),pch=15,bty="n",cex=1)
			par(xpd=FALSE)
			nb.W <- length(unique(object@zc[[id]]))

			sum.rho <- 0
			if(nb.W!=1){
				for(i in 1:(nb.W-1)){
					rho <- length(which(object@zc[[id]]==i))/length(object@zc[[id]])
					sum.rho <- sum.rho + rho
					abline(v=sum.rho,lwd=3, col="red")
				}
			}
			

			sum.gammas <- 0
			for(i in 1:(nb.V-1)){
				sum.gammas <- sum.gammas + gammas[i]
				abline(h=sum.gammas,lwd=3, col="red")
			}

		}
	}
	if(object@name == "Clust"){
		D <- length(m)
		par(mfrow=c(1,D))
		for(id in 1:D){
			tmp <- object@xhat[[id]][sort(object@zr,index.return=TRUE)$ix,1:ncol(object@xhat[[id]])]
			image(t(tmp),xaxt='n',yaxt='n', main='clustering',cex.main=1.5,col = gray(m[id]:1/m[id]))
		    par(xpd=TRUE, mar=c(2,1,4,7))

			legend(x=1.01,y=0.65,legend = 1:m[id],col = gray(m[id]:1/m[id]),pch=15,bty="n",cex=1)
			par(xpd=FALSE)
			sum.gammas <- 0
			for(i in 1:(nb.V-1)){
				sum.gammas <- sum.gammas + gammas[i]
				abline(h=sum.gammas,lwd=3, col="red")
			}
		}
	}
	
}

bosplot <- function(object){
	par(xpd=TRUE, mar=c(2,1,4,7))
	par(xpd=FALSE)

	m <- object@m

	nb.V <- length(unique(object@zr))
	gammas <- rep(0,nb.V)
	for(i in 1:nb.V){
		gammas[i] <- length(which(object@zr==i))/length(object@zr)
	}


	
	if(object@name == "ClassifM"){
		D <- length(m)
		par(mfrow=c(1,D))
		for(id in 1:D){
			tmp <- object@xhat[[id]][sort(object@zr,index.return=TRUE)$ix,1:ncol(object@xhat[[id]])]
			par(xpd=TRUE, mar=c(2,1,4,7))
			image(t(tmp),xaxt='n',yaxt='n', main='classification',cex.main=1.5,col = gray(m[id]:1/m[id]))


			legend(x=1.01,y=0.65,legend = 1:m[id],col = gray(m[id]:1/m[id]),pch=15,bty="n",cex=1)
			par(xpd=FALSE)
			sum.gammas <- 0
			for(i in 1:(nb.V-1)){
				sum.gammas <- sum.gammas + gammas[i]
				abline(h=sum.gammas,lwd=3, col="red")
			}
		}
	}
	if(object@name == "Classif"){
		D <- length(m)
		par(mfrow=c(1,D))
		for(id in 1:D){
			tmp <- object@xhat[[id]][sort(object@zr,index.return=TRUE)$ix,
					sort(object@zc[[id]],index.return=TRUE)$ix]
			par(xpd=TRUE, mar=c(2,1,4,7))
			image(t(tmp),xaxt='n',yaxt='n', main='classification',cex.main=1.5,col = gray(m[id]:1/m[id]))


			legend(x=1.01,y=0.65,legend = 1:m[id],col = gray(m[id]:1/m[id]),pch=15,bty="n",cex=1)
			par(xpd=FALSE)
			nb.W <- length(unique(object@zc[[id]]))

			sum.rho <- 0
			if(nb.W!=1){
				for(i in 1:(nb.W-1)){
					rho <- length(which(object@zc[[id]]==i))/length(object@zc[[id]])
					sum.rho <- sum.rho + rho
					abline(v=sum.rho,lwd=3, col="red")
				}
			}
			

			sum.gammas <- 0
			for(i in 1:(nb.V-1)){
				sum.gammas <- sum.gammas + gammas[i]
				abline(h=sum.gammas,lwd=3, col="red")
			}
		}
	}
	if(object@name == "Coclust"){
		D <- length(m)
		par(mfrow=c(1,D))
		for(id in 1:D){
			tmp <- object@xhat[[id]][sort(object@zr,index.return=TRUE)$ix,
					sort(object@zc[[id]],index.return=TRUE)$ix]
			par(xpd=TRUE, mar=c(2,1,4,7))
			image(t(tmp),xaxt='n',yaxt='n', main='co-clustering',cex.main=1.5,col = gray(m[id]:1/m[id]))


			legend(x=1.01,y=0.65,legend = 1:m[id],col = gray(m[id]:1/m[id]),pch=15,bty="n",cex=1)
			par(xpd=FALSE)
			nb.W <- length(unique(object@zc[[id]]))

			sum.rho <- 0
			if(nb.W!=1){
				for(i in 1:(nb.W-1)){
					rho <- length(which(object@zc[[id]]==i))/length(object@zc[[id]])
					sum.rho <- sum.rho + rho
					abline(v=sum.rho,lwd=3, col="red")
				}
			}
			

			sum.gammas <- 0
			for(i in 1:(nb.V-1)){
				sum.gammas <- sum.gammas + gammas[i]
				abline(h=sum.gammas,lwd=3, col="red")
			}

		}
	}
	if(object@name == "Clust"){
		D <- length(m)
		par(mfrow=c(1,D))
		for(id in 1:D){
			tmp <- object@xhat[[id]][sort(object@zr,index.return=TRUE)$ix,1:ncol(object@xhat[[id]])]
			image(t(tmp),xaxt='n',yaxt='n', main='clustering',cex.main=1.5,col = gray(m[id]:1/m[id]))
		    par(xpd=TRUE, mar=c(2,1,4,7))

			legend(x=1.01,y=0.65,legend = 1:m[id],col = gray(m[id]:1/m[id]),pch=15,bty="n",cex=1)
			par(xpd=FALSE)
			sum.gammas <- 0
			for(i in 1:(nb.V-1)){
				sum.gammas <- sum.gammas + gammas[i]
				abline(h=sum.gammas,lwd=3, col="red")
			}
		}
	}
	
}