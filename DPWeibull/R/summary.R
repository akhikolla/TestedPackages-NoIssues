summary.dpm<-function(object,...){

    ans<-NULL
    if(object$usertime){
    index<-floor(c(1/4,1/2,3/4,1)*length(object$Spred))
    ans<-rbind(object$predtime[index], object$Spred[index],object$Spredu[index],object$Spredl[index])
    }else{
    index<-c(11,21,31,41)
    ans<-rbind(object$predtime[index], object$Spred[index],object$Spredu[index],object$Spredl[index])
    }
    rownames(ans)<-c("Time","Estimated Survival","Upper Credible Interval", "Lower Credible Interval")
    colnames(ans)<-c("1/4 Observation Time","1/2 Observation Time","3/4 Observation Time","Observation Time")
    class(ans)<-"summary.dpm"
    ans
}

summary.ddp<-function(object,...){
    ans<-NULL
    if(object$usertime){
    	    index<-floor(c(1/4,1/2,3/4,1)*ncol(object$loghr.est))
    }else{
    	    index<-c(11,21,31,41)
    }
	    for(i in 1:nrow(object$loghr.est)){
	    ans[[i]]<-rbind(object$predtime[index],object$loghr.est[i,index], object$loghru[i,index],object$loghrl[i,index])
	    rownames(ans[[i]])<-c("Time","Estimated Log Hazard Ratio","Upper Credible Interval", "Lower Credible Interval")
	    colnames(ans[[i]])<-c("1/4 Observation Time","1/2 Observation Time","3/4 Observation Time","Observation Time")
	    }
   names(ans)<-object$covnames
    class(ans)<-"summary.ddp"
    ans
}

summary.dpmcomp<-function(object,...){
    ans<-NULL
    if(object$usertime){
    index<-floor(c(1/4,1/2,3/4,1)*length(object$CIF1.est))
    
    }else{
    index<-c(11,21,31,41)
    
    }
    ans$CIF1<-rbind(object$predtime[index], object$CIF1.est[index],object$CIF1u[index],object$CIF1l[index])
    ans$CIF2<-rbind(object$predtime[index], object$CIF2.est[index],object$CIF2u[index],object$CIF2l[index])
    rownames(ans$CIF1)<-c("Time","Estimated Cumulative Incidence Function for Cause 1","Upper Credible Interval", "Lower Credible Interval")

    colnames(ans$CIF1)<-c("1/4 Observation Time","1/2 Observation Time","3/4 Observation Time","Observation Time")

    rownames(ans$CIF2)<-c("Time","Estimated Cumulative Incidence Function for Cause 2","Upper Credible Interval", "Lower Credible Interval")

    colnames(ans$CIF2)<-c("1/4 Observation Time","1/2 Observation Time","3/4 Observation Time","Observation Time")
    class(ans)<-"summary.dpmcomp"
    ans
}

summary.ddpcomp<-function(object,...){
    ans<-NULL
    if(object$usertime){
    	    index<-floor(c(1/4,1/2,3/4,1)*ncol(object$loghr.est))
    }else{
    	    index<-c(11,21,31,41)
	}
	    for(i in 1:nrow(object$loghr.est)){
	    ans[[i]]<-rbind(object$predtime[index],object$loghr.est[i,index], object$loghru[i,index],object$loghrl[i,index])
	    rownames(ans[[i]])<-c("Time","Estimated Log Hazard Ratio","Upper Credible Interval", "Lower Credible Interval")
	    colnames(ans[[i]])<-c("1/4 Observation Time","1/2 Observation Time","3/4 Observation Time","Observation Time")
	    }
 	    
   names(ans)<-object$covnames
    class(ans)<-"summary.ddpcomp"
    ans
}

summary.predddp<-function(object,...){
    ans<-NULL
    for(i in 1:nrow(object$Spred)){
    index<-floor(c(1/4,1/2,3/4,1)*ncol(object$Spred))
    ans$S[[i]]<-rbind(object$tpred[index], object$Spred[i,index],object$Spredu[i,index],object$Spredl[i,index])
    rownames(ans$S[[i]])<-c("Time","Estimated Survival","Upper Credible Interval", "Lower Credible Interval")
    colnames(ans$S[[i]])<-c("1/4 Prediction Time","1/2 Prediction Time","3/4 Prediction Time","Prediction Time")
    }
    class(ans)<-"summary.predddp"
    ans
}


summary.predddpcomp<-function(object,...){
    ans<-NULL
    for(i in 1:nrow(object$Fpred)){
    index<-floor(c(1/4,1/2,3/4,1)*ncol(object$Fpred))
    ans$CIF1[[i]]<-rbind(object$tpred[index], object$Fpred[index],object$Fpredu[index],object$Fpredl[index])
    rownames(ans$CIF1[[i]])<-c("Time","Estimated Cumulative Incidence Function","Upper Credible Interval", "Lower Credible Interval")
    colnames(ans$CIF1[[i]])<-c("1/4 Prediction Time","1/2 Prediction Time","3/4 Prediction Time","Prediction Time")
    }
    class(ans)<-"summary.predddpcomp"
    ans
}
