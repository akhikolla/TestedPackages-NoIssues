load("data/survmedian.rda")
survmedian<-survmedian[order(survmedian)]

findhighpct<-function(p){

whichint<-100-findInterval(p,survmedian)
100/whichint
}

is.binary <- function(v) {
  x <- unique(v)
  length(x) - sum(is.na(x)) == 2L
}


dpmdiffalpha<-function(alpha,result){
   result$Spredu<-apply(result$S,2,quantile,1-alpha/2)
   result$Spredl<-apply(result$S,2,quantile,alpha/2)

   result$dpredu<-apply(result$d,2,quantile,1-alpha/2)
   result$dpredl<-apply(result$d,2,quantile,alpha/2)

   result$hpredu<-apply(result$h,2,quantile,1-alpha/2)
   result$hpredl<-apply(result$h,2,quantile,alpha/2)
   result
}

ddpcomppreddiffalpha<-function(alpha,result){
   result$Spredu<-apply(result$S,2,quantile,1-alpha/2)
   result$Spredl<-apply(result$S,2,quantile,alpha/2)

   result$dpredu<-apply(result$d,2,quantile,1-alpha/2)
   result$dpredl<-apply(result$d,2,quantile,alpha/2)

   result$hpredu<-apply(result$h,2,quantile,1-alpha/2)
   result$hpredl<-apply(result$h,2,quantile,alpha/2)
   result
}

ddpdiffalpha<-function(alpha,result){
  result$loghrl<-matrix(apply(result$loghr,2,quantile,alpha/2),byrow=TRUE,nrow=ncol(result$x))/result$xscale
  result$loghru<-matrix(apply(result$loghr,2,quantile,1-alpha/2),byrow=TRUE,nrow=ncol(result$x))/result$xscale
   result
}

dpmcompdiffalpha<-function(alpha,result){
   result$CIF1u<-apply(result$CIF1,2,quantile,1-alpha/2)
   result$CIF1l<-apply(result$CIF1,2,quantile,alpha/2)

   result$d1u<-apply(result$d1,2,quantile,1-alpha/2)
   result$d1l<-apply(result$d1,2,quantile,alpha/2)

   result$h1u<-apply(result$h1,2,quantile,1-alpha/2)
   result$h1l<-apply(result$h1,2,quantile,alpha/2)

   result$CIF2u<-apply(result$CIF2,2,quantile,1-alpha/2)
   result$CIF2l<-apply(result$CIF2,2,quantile,alpha/2)

   result$d2u<-apply(result$d2,2,quantile,1-alpha/2)
   result$d2l<-apply(result$d2,2,quantile,alpha/2)

   result$h2u<-apply(result$h2,2,quantile,1-alpha/2)
   result$h2l<-apply(result$h2,2,quantile,alpha/2)
   result
}

confband<-function(alpha,result){
	resultcopy<-result
	for(i in 1:ncol(result)){
		resultcopy[,i]<-result[order(-result[,i]),i]
	}
	minindex<-ceiling(nrow(result)*alpha/2)
	perc<-1:minindex
	for(i in 1:minindex){
                 temp<-0
                 for(j in 1:nrow(result)){
		temp<-temp+as.numeric(any(result[j,]>resultcopy[i,])|any(result[j,]<resultcopy[1+nrow(result)-i,]))
		}
		perc[i]<-temp/nrow(result)
	}
       percindex<-max(which(perc<alpha))

       rbind(resultcopy[percindex,],resultcopy[1+nrow(result)-percindex,])
 }
