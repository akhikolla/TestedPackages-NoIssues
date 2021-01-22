
plot.dpm<-function(x,simultaneous=FALSE,...){
	if((x$simultaneous==FALSE)&(simultaneous==FALSE)){
  	plot(c(0,x$predtime),c(1,x$Spred),type="l",lwd=3,main="Survival",ylab="",xlab="Time",ylim=c(min(x$Spredl,na.rm=TRUE),1))
	}else{
	  x$Sbandl<-confband(x$alpha,x$S)[1,]
	  x$Sbandu<-confband(x$alpha,x$S)[2,]
  	plot(c(0,x$predtime),c(1,x$Spred),type="l",lwd=3,main="Survival",ylab="",xlab="Time",ylim=c(min(x$Sbandl,na.rm=TRUE),1))
  	lines(c(0,x$predtime),c(1,x$Sbandu),lty=3,lwd=3)
  	lines(c(0,x$predtime),c(1,x$Sbandl),lty=3,lwd=3)
	}
  lines(c(0,x$predtime),c(1,x$Spredu),lty=2,lwd=3)
  lines(c(0,x$predtime),c(1,x$Spredl),lty=2,lwd=3)
  
	if((x$simultaneous==FALSE)&(simultaneous==FALSE)){
  	plot(c(0,x$predtime),c(0,x$dpred),type="l",lwd=3,main="Density",ylab="",xlab="Time",ylim=c(0,max(x$dpredu,na.rm=TRUE)))
	}else{
  	x$dbandl<-confband(x$alpha,x$d)[1,]
  	x$dbandu<-confband(x$alpha,x$d)[2,]
  	plot(c(0,x$predtime),c(0,x$dpred),type="l",lwd=3,main="Density",ylab="",xlab="Time",ylim=c(0,max(x$dbandu,na.rm=TRUE)))
  	lines(c(0,x$predtime),c(0,x$dbandu),lty=3,lwd=3)
  	lines(c(0,x$predtime),c(0,x$dbandl),lty=3,lwd=3)
	}
  lines(c(0,x$predtime),c(0,x$dpredu),lty=2,lwd=3)
  lines(c(0,x$predtime),c(0,x$dpredl),lty=2,lwd=3)
    
	if((x$simultaneous==FALSE)&(simultaneous==FALSE)){
  	plot(x$predtime,x$hpred,type="l",lwd=3,main="Hazard",ylab="",xlab="Time",ylim=c(min(x$hpredl,na.rm=TRUE),max(x$hpredu,na.rm=TRUE)))
  	lines(x$predtime,x$hpredu,lty=2,lwd=3)
  	lines(x$predtime,x$hpredl,lty=2,lwd=3)
	}else{
	  x$hbandl<-confband(x$alpha,x$h)[1,]
	  x$hbandu<-confband(x$alpha,x$h)[2,]
  	plot(x$predtime,x$hpred,type="l",lwd=3,main="Hazard",ylab="",xlab="Time",ylim=c(min(x$hbandl,na.rm=TRUE),max(x$hbandu,na.rm=TRUE)))
  	lines(x$predtime,x$hbandu,lty=3,lwd=3)
  	lines(x$predtime,x$hbandl,lty=3,lwd=3)
	}
  lines(x$predtime,x$hpredu,lty=2,lwd=3)
  lines(x$predtime,x$hpredl,lty=2,lwd=3)
}

plot.dpmcomp<-function(x,simultaneous=FALSE,...){
	if((x$simultaneous==FALSE)&(simultaneous==FALSE)){
  	plot(c(0,x$predtime),c(0,x$CIF1.est),type="l",col="red",lwd=3,main="Cumulative Incidence Functions",ylab="",xlab="Time",
  	     ylim=c(0,max(c(x$CIF1u,x$CIF2u),na.rm=TRUE)))
  	lines(c(0,x$predtime),c(0,x$CIF2.est),lwd=3,col="blue")
  	legend("topleft",c("Event 1", "Event 2"), lwd=c(3,3), lty=c(1,1), col=c("red", "blue"))
	}else{
	  x$CIF1bandl<-confband(x$alpha,x$CIF1)[1,]
	  x$CIF1bandu<-confband(x$alpha,x$CIF1)[2,]
	  x$CIF2bandl<-confband(x$alpha,x$CIF2)[1,]
	  x$CIF2bandu<-confband(x$alpha,x$CIF2)[2,]
  	plot(c(0,x$predtime),c(0,x$CIF1.est),type="l",col="red",lwd=3,main="Cumulative Incidence Functions",ylab="",xlab="Time",
  	     ylim=c(0,max(c(x$CIF1bandu,x$CIF2bandu),na.rm=TRUE)))
  	lines(c(0,x$predtime),c(0,x$CIF1bandu),lty=3,lwd=3,col="red")
  	lines(c(0,x$predtime),c(0,x$CIF1bandl),lty=3,lwd=3,col="red")
  	lines(c(0,x$predtime),c(0,x$CIF2bandu),lty=3,lwd=3,col="blue")
  	lines(c(0,x$predtime),c(0,x$CIF2bandl),lty=3,lwd=3,col="blue")
	  legend("bottomright",c("Event 1", "Event 2"), lwd=c(3,3), lty=c(1,1), col=c("red", "blue"))
	}
  lines(c(0,x$predtime),c(0,x$CIF2.est),lwd=3,col="blue")
  lines(c(0,x$predtime),c(0,x$CIF1u),lty=2,lwd=3,col="red")
  lines(c(0,x$predtime),c(0,x$CIF1l),lty=2,lwd=3,col="red")
  lines(c(0,x$predtime),c(0,x$CIF2u),lty=2,lwd=3,col="blue")
  lines(c(0,x$predtime),c(0,x$CIF2l),lty=2,lwd=3,col="blue")
  
	if((x$simultaneous==FALSE)&(simultaneous==FALSE)){
  	plot(c(0,x$predtime),c(0,x$d1.est),type="l",col="red",lwd=3,main="Subdistribution Density Functions",ylab="",xlab="Time",
  	     ylim=c(0,max(c(x$d1u,x$d2u),na.rm=TRUE)))
	}else{
	  x$d1bandl<-confband(x$alpha,x$d1)[1,]
	  x$d1bandu<-confband(x$alpha,x$d1)[2,]
	  x$d2bandl<-confband(x$alpha,x$d2)[1,]
	  x$d2bandu<-confband(x$alpha,x$d2)[2,]
  	plot(c(0,x$predtime),c(0,x$d1.est),type="l",col="red",lwd=3,main="Subdistribution Density Functions",ylab="",xlab="Time",
  	     ylim=c(0,max(c(x$d1bandu,x$d2bandu),na.rm=TRUE)))
  	lines(c(0,x$predtime),c(0,x$d1bandu),lty=3,lwd=3,col="red")
  	lines(c(0,x$predtime),c(0,x$d1bandl),lty=3,lwd=3,col="red")
  	lines(c(0,x$predtime),c(0,x$d2bandu),lty=3,lwd=3,col="blue")
  	lines(c(0,x$predtime),c(0,x$d2bandl),lty=3,lwd=3,col="blue")
	}
  lines(c(0,x$predtime),c(0,x$d2.est),lwd=3,col="blue")
  lines(c(0,x$predtime),c(0,x$d1u),lty=2,lwd=3,col="red")
  lines(c(0,x$predtime),c(0,x$d1l),lty=2,lwd=3,col="red")
  lines(c(0,x$predtime),c(0,x$d2u),lty=2,lwd=3,col="blue")
  lines(c(0,x$predtime),c(0,x$d2l),lty=2,lwd=3,col="blue")
  
	if((x$simultaneous==FALSE)&(simultaneous==FALSE)){
  	plot(x$predtime,x$h1.est,type="l",col="red",lwd=3,main="Subdistribution Hazard Functions",ylab="",xlab="Time",
  	     ylim=c(min(c(x$h1l,x$h2l),na.rm=TRUE),max(c(x$h1u,x$h2u),na.rm=TRUE)))
	}else{
	  x$h1bandl<-confband(x$alpha,x$h1)[1,]
	  x$h1bandu<-confband(x$alpha,x$h1)[2,]
	  x$h2bandl<-confband(x$alpha,x$h2)[1,]
	  x$h2bandu<-confband(x$alpha,x$h2)[2,]
  	plot(x$predtime,x$h1.est,type="l",col="red",lwd=3,main="Subdistribution Hazard Functions",ylab="",xlab="Time",
  	     ylim=c(min(c(x$h1bandl,x$h2bandl),na.rm=TRUE),max(c(x$h1bandu,x$h2bandu),na.rm=TRUE)))
  	lines(x$predtime,x$h1bandu,lty=3,lwd=3,col="red")
  	lines(x$predtime,x$h1bandl,lty=3,lwd=3,col="red")
  	lines(x$predtime,x$h2bandu,lty=3,lwd=3,col="blue")
  	lines(x$predtime,x$h2bandl,lty=3,lwd=3,col="blue")
	}
  lines(x$predtime,x$h1u,lty=2,lwd=3,col="red")
  lines(x$predtime,x$h1l,lty=2,lwd=3,col="red")
  lines(x$predtime,x$h2.est,lwd=3,col="blue")
  lines(x$predtime,x$h2u,lty=2,lwd=3,col="blue")
  lines(x$predtime,x$h2l,lty=2,lwd=3,col="blue")
}

plot.ddp<-function(x,simultaneous=FALSE,exp=FALSE,...){
  if((x$simultaneous==FALSE)&(simultaneous==TRUE)){
  	x$loghrbandl<-matrix(confband(x$alpha,x$loghr)[1,],byrow=TRUE,nrow=ncol(x$x))/x$xscale
    x$loghrbandu<-matrix(confband(x$alpha,x$loghr)[2,],byrow=TRUE,nrow=ncol(x$x))/x$xscale
  	if(exp==TRUE){
  		x$hrbandl<-exp(x$loghrbandl)
  		x$hrbandu<-exp(x$loghrbandu)
  	}
	}
	if((x$simultaneous==FALSE)&(simultaneous==FALSE)){
		if(exp==FALSE){
  		for(i in 1:nrow(x$loghr.est)){
    		ytop<-max(x$loghru,na.rm=TRUE)
    		ybot<-min(x$loghrl,na.rm=TRUE)
    		padName<-x$covnames[i]
    		if(grepl("factor",padName,fixed=TRUE)){
    		  factorName<-sub(").*", ")", padName)
    		  # now called factor(Stage)
    		  reference<-unlist(x$xlevels[factorName])[1]
    		  padName<-gsub("factor[(]", "", padName)
    		  padName<-gsub("[)]", "=", padName)
    		  factorName<-gsub("factor[(]", "", factorName)
    		  factorName<-gsub("[)]", "=", factorName)
    		  reference<-paste(factorName,reference)
    		  padName<-paste(padName,"vs",reference)
    		}
    		plot(x$predtime,x$loghr.est[i,],type="l",lwd=3,main=paste("Log Hazard Ratio over Time for Covariate ",padName,sep="")
    		 ,ylab="",xlab="Time",ylim=c(ybot,ytop))
    		lines(x$predtime,x$loghrl[i,],lty=2,lwd=3)
    		lines(x$predtime,x$loghru[i,],lty=2,lwd=3)
  		}
		}else{
      x$hr.est<-exp(x$loghr.est)
      x$hrl<-exp(x$loghrl)
      x$hru<-exp(x$loghru)
    	for(i in 1:nrow(x$loghr.est)){
    		ytop<-max(x$hru,na.rm=TRUE)
    		ybot<-min(x$hrl,na.rm=TRUE)
    		padName<-x$covnames[i]
    		if(grepl("factor",padName,fixed=TRUE)){
    		  factorName<-sub(").*", ")", padName)
    		  # now called factor(Stage)
    		  reference<-unlist(x$xlevels[factorName])[1]
    		  padName<-gsub("factor[(]", "", padName)
    		  padName<-gsub("[)]", "=", padName)
    		  factorName<-gsub("factor[(]", "", factorName)
    		  factorName<-gsub("[)]", "=", factorName)
    		  reference<-paste(factorName,reference)
    		  padName<-paste(padName,"vs",reference)
    		}
    		plot(x$predtime,x$hr.est[i,],type="l",lwd=3,
    		main=paste("Hazard Ratio over Time for Covariate ",padName,sep="")
    		 ,ylab="",xlab="Time",ylim=c(ybot,ytop))
    		lines(x$predtime,x$hrl[i,],lty=2,lwd=3)
    		lines(x$predtime,x$hru[i,],lty=2,lwd=3)
		  }
		}
	}else{
		if(exp==FALSE){
  		for(i in 1:nrow(x$loghr.est)){
    		ytop<-max(x$loghr,na.rm=TRUE)
    		ybot<-min(x$loghr,na.rm=TRUE)
    		padName<-x$covnames[i]
    		if(grepl("factor",padName,fixed=TRUE)){
    		  factorName<-sub(").*", ")", padName)
    		  # now called factor(Stage)
    		  reference<-unlist(x$xlevels[factorName])[1]
    		  padName<-gsub("factor[(]", "", padName)
    		  padName<-gsub("[)]", "=", padName)
    		  factorName<-gsub("factor[(]", "", factorName)
    		  factorName<-gsub("[)]", "=", factorName)
    		  reference<-paste(factorName,reference)
    		  padName<-paste(padName,"vs",reference)
    		}
    		plot(x$predtime,x$loghr.est[i,],type="l",lwd=3,main=paste("Log Hazard Ratio over Time for Covariate ",padName,sep="")
    		 ,ylab="",xlab="Time",ylim=c(ybot,ytop))
    		lines(x$predtime,x$loghrl[i,],lty=2,lwd=3)
    		lines(x$predtime,x$loghru[i,],lty=2,lwd=3)
    		lines(x$predtime,x$loghrbandl[i,],lty=3,lwd=3)
    		lines(x$predtime,x$loghrbandu[i,],lty=3,lwd=3)
  		}
		}else{
      x$hr<-exp(x$loghr)
      x$hr.est<-exp(x$loghr.est)
      x$hrl<-exp(x$loghrl)
      x$hru<-exp(x$loghru)
  		for(i in 1:nrow(x$loghr.est)){
    		ytop<-max(x$hr,na.rm=TRUE)
    		ybot<-min(x$hr,na.rm=TRUE)
    		padName<-x$covnames[i]
    		if(grepl("factor",padName,fixed=TRUE)){
    		  factorName<-sub(").*", ")", padName)
    		  # now called factor(Stage)
    		  reference<-unlist(x$xlevels[factorName])[1]
    		  padName<-gsub("factor[(]", "", padName)
    		  padName<-gsub("[)]", "=", padName)
    		  factorName<-gsub("factor[(]", "", factorName)
    		  factorName<-gsub("[)]", "=", factorName)
    		  reference<-paste(factorName,reference)
    		  padName<-paste(padName,"vs",reference)
    		}
    		plot(x$predtime,x$hr.est[i,],type="l",lwd=3,
    		main=paste("Hazard Ratio over Time for Covariate ",padName,sep="")
    		 ,ylab="",xlab="Time",ylim=c(ybot,ytop))
    		lines(x$predtime,x$hrl[i,],lty=2,lwd=3)
    		lines(x$predtime,x$hru[i,],lty=2,lwd=3)
    		lines(x$predtime,x$hrbandl[i,],lty=3,lwd=3)
    		lines(x$predtime,x$hrbandu[i,],lty=3,lwd=3)
  		}
		}
	}
}

plot.ddpcomp<-function(x,simultaneous=FALSE,exp=FALSE,...){
  if((x$simultaneous==FALSE)&(simultaneous==TRUE)){
	  x$loghrbandl<-matrix(confband(x$alpha,x$loghr)[1,],byrow=TRUE,nrow=ncol(x$x))/x$xscale
  	x$loghrbandu<-matrix(confband(x$alpha,x$loghr)[2,],byrow=TRUE,nrow=ncol(x$x))/x$xscale
	  if(exp==TRUE){
  		x$hrbandl<-exp(x$loghrbandl)
  		x$hrbandu<-exp(x$loghrbandu)
	  }
  }
	if((x$simultaneous==FALSE)&(simultaneous==FALSE)){
		if(exp==FALSE){
  		for(i in 1:nrow(x$loghr.est)){
    		ytop<-max(x$loghru,na.rm=TRUE)
    		ybot<-min(x$loghrl,na.rm=TRUE)
    		padName<-x$covnames[i]
    		if(grepl("factor",padName,fixed=TRUE)){
    		  factorName<-sub(").*", ")", padName)
    		  # now called factor(Stage)
    		  reference<-unlist(x$xlevels[factorName])[1]
    		  padName<-gsub("factor[(]", "", padName)
    		  padName<-gsub("[)]", "=", padName)
    		  factorName<-gsub("factor[(]", "", factorName)
    		  factorName<-gsub("[)]", "=", factorName)
    		  reference<-paste(factorName,reference)
    		  padName<-paste(padName,"vs",reference)
    		}
    		plot(x$predtime,x$loghr.est[i,],type="l",lwd=3,
    		main=paste("Log Subdistribution Hazard Ratio of \n Covariate ",padName," for Event 1",sep="")
    		 ,ylab="",xlab="Time",ylim=c(ybot,ytop))
    		lines(x$predtime,x$loghrl[i,],lty=2,lwd=3)
    		lines(x$predtime,x$loghru[i,],lty=2,lwd=3)
	  	}
	  }else{
      x$hr.est<-exp(x$loghr.est)
      x$hrl<-exp(x$loghrl)
      x$hru<-exp(x$loghru)
  		for(i in 1:nrow(x$loghr.est)){
    		ytop<-max(x$hru,na.rm=TRUE)
    		ybot<-min(x$hrl,na.rm=TRUE)
    		padName<-x$covnames[i]
    		if(grepl("factor",padName,fixed=TRUE)){
    		  factorName<-sub(").*", ")", padName)
    		  # now called factor(Stage)
    		  reference<-unlist(x$xlevels[factorName])[1]
    		  padName<-gsub("factor[(]", "", padName)
    		  padName<-gsub("[)]", "=", padName)
    		  factorName<-gsub("factor[(]", "", factorName)
    		  factorName<-gsub("[)]", "=", factorName)
    		  reference<-paste(factorName,reference)
    		  padName<-paste(padName,"vs",reference)
    		}
    		plot(x$predtime,x$hr.est[i,],type="l",lwd=3,
    		main=paste("Subdistribution Hazard Ratio of \n Covariate ",padName," for Event 1",sep="")
    		 ,ylab="",xlab="Time",ylim=c(ybot,ytop))
    		lines(x$predtime,x$hrl[i,],lty=2,lwd=3)
    		lines(x$predtime,x$hru[i,],lty=2,lwd=3)
  		}
		}
  }else{
		if(exp==FALSE){
  		for(i in 1:nrow(x$loghr.est)){
    		ytop<-max(x$loghr,na.rm=TRUE)
    		ybot<-min(x$loghr,na.rm=TRUE)
    		padName<-x$covnames[i]
    		if(grepl("factor",padName,fixed=TRUE)){
    		  factorName<-sub(").*", ")", padName)
    		  # now called factor(Stage)
    		  reference<-unlist(x$xlevels[factorName])[1]
    		  padName<-gsub("factor[(]", "", padName)
    		  padName<-gsub("[)]", "=", padName)
    		  factorName<-gsub("factor[(]", "", factorName)
    		  factorName<-gsub("[)]", "=", factorName)
    		  reference<-paste(factorName,reference)
    		  padName<-paste(padName,"vs",reference)
    		}
    		plot(x$predtime,x$loghr.est[i,],type="l",lwd=3,
    		main=paste("Log Subdistribution Hazard Ratio of \n Covariate ",padName," for Event 1",sep="")
    		 ,ylab="",xlab="Time",ylim=c(ybot,ytop))
    		lines(x$predtime,x$loghrl[i,],lty=2,lwd=3)
    		lines(x$predtime,x$loghru[i,],lty=2,lwd=3)
    		lines(x$predtime,x$loghrbandl[i,],lty=3,lwd=3)
    		lines(x$predtime,x$loghrbandu[i,],lty=3,lwd=3)
  		}
		}else{
      x$hr<-exp(x$loghr)
      x$hr.est<-exp(x$loghr.est)
      x$hrl<-exp(x$loghrl)
      x$hru<-exp(x$loghru)
  		for(i in 1:nrow(x$loghr.est)){
    		ytop<-max(x$hr,na.rm=TRUE)
    		ybot<-min(x$hr,na.rm=TRUE)
    		padName<-x$covnames[i]
    		if(grepl("factor",padName,fixed=TRUE)){
    		  factorName<-sub(").*", ")", padName)
    		  # now called factor(Stage)
    		  reference<-unlist(x$xlevels[factorName])[1]
    		  padName<-gsub("factor[(]", "", padName)
    		  padName<-gsub("[)]", "=", padName)
    		  factorName<-gsub("factor[(]", "", factorName)
    		  factorName<-gsub("[)]", "=", factorName)
    		  reference<-paste(factorName,reference)
    		  padName<-paste(padName,"vs",reference)
    		}
    		plot(x$predtime,x$hr.est[i,],type="l",lwd=3,
    		main=paste("Log Subdistribution Hazard Ratio of \n Covariate ",padName," for Event 1",sep="")
    		 ,ylab="",xlab="Time",ylim=c(ybot,ytop))
    		lines(x$predtime,x$hrl[i,],lty=2,lwd=3)
    		lines(x$predtime,x$hru[i,],lty=2,lwd=3)
    		lines(x$predtime,x$hrbandl[i,],lty=3,lwd=3)
    		lines(x$predtime,x$hrbandu[i,],lty=3,lwd=3)
  		}
		}
	}
}

plot.predddpcomp<-function(x,...){
	for(i in 1:nrow(x$Fpred)){
	  plot(c(0,x$tpred),c(0,x$Fpred[i,]),
  	main=paste("Cumulative Incidence Function Estimate\n with New Data ",i, " for Event 1", sep=""),
  	type="l",lwd=3,xlab="Time",ylab="",ylim=c(min(x$Fpredl,na.rm=TRUE),max(x$Fpredu,na.rm=TRUE)))
	  lines(c(0,x$tpred),c(0,x$Fpredl[i,]),lwd=3,lty=2)
	  lines(c(0,x$tpred),c(0,x$Fpredu[i,]),lwd=3,lty=2)
	  
	  plot(c(0,x$tpred),c(0,x$dpred[i,]),
  	main=paste("Cause-specific Density Estimate\n with New Data ",i," for Event 1", sep=""),
  	type="l",lwd=3,xlab="Time",ylab="",ylim=c(min(x$dpredl,na.rm=TRUE),max(x$dpredu,na.rm=TRUE)))
	  lines(c(0,x$tpred),c(0,x$dpredl[i,]),lwd=3,lty=2)
	  lines(c(0,x$tpred),c(0,x$dpredu[i,]),lwd=3,lty=2)
	  
	  plot(x$tpred,x$hpred[i,],
  	main=paste("Subdistribution Hazard Estimate\n with New Data ",i," for Event 1", sep=""),
  	type="l",lwd=3,xlab="Time",ylab="",ylim=c(min(x$hpredl,na.rm=TRUE),max(x$hpredu,na.rm=TRUE)))
	  lines(x$tpred,x$hpredl[i,],lwd=3,lty=2)
	  lines(x$tpred,x$hpredu[i,],lwd=3,lty=2)
	}
}



plot.predddp<-function(x,...){
  for(i in 1:nrow(x$Spred)){
      plot(c(0,x$tpred),c(1,x$Spred[i,]),
      main=paste("Survival Estimate with New Data ",i,  sep=""),
      type="l",lwd=3,xlab="Time",ylab="",ylim=c(min(x$Spredl,na.rm=TRUE),max(x$Spredu,na.rm=TRUE)))
      lines(c(0,x$tpred),c(1,x$Spredl[i,]),lwd=3,lty=2)
      lines(c(0,x$tpred),c(1,x$Spredu[i,]),lwd=3,lty=2)
      
      plot(c(0,x$tpred),c(0,x$dpred[i,]),
      main=paste("Density Estimate with New Data ",i, sep=""),
      type="l",lwd=3,xlab="Time",ylab="",ylim=c(min(x$dpredl,na.rm=TRUE),max(x$dpredu,na.rm=TRUE)))
      lines(c(0,x$tpred),c(0,x$dpredl[i,]),lwd=3,lty=2)
      lines(c(0,x$tpred),c(0,x$dpredu[i,]),lwd=3,lty=2)
      
      plot(x$tpred,x$hpred[i,],
      main=paste("Hazard Estimate with New Data ",i,sep=""),
      type="l",lwd=3,xlab="Time",ylab="",ylim=c(min(x$hpredl,na.rm=TRUE),max(x$hpredu,na.rm=TRUE)))
      lines(x$tpred,x$hpredl[i,],lwd=3,lty=2)
      lines(x$tpred,x$hpredu[i,],lwd=3,lty=2)
  }
}


