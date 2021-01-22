predict.ddp<-function(object,newdata,alpha=0.05,tpred=NULL,...){
if(is.null(tpred)){
tpred<-object$predtime
}
tpred<-tpred/object$high.pct*10

xpred<-model.matrix(object$Terms,newdata,xlev = object$xlevels)
xpred<-xpred[,-1]
xpred<-(xpred-matrix(rep(object$xmean,times=nrow(xpred)),nrow=nrow(xpred), byrow=TRUE))/matrix(rep(2*object$xsd,times=nrow(xpred)),nrow=nrow(xpred), byrow=TRUE)	
predobject<-.Call('DPWeibull_predreg', PACKAGE = 'DPWeibull',
object$alpharec,object$lambdascaled,object$betarec,xpred,tpred,alpha)
class(predobject)<-"predddp"
predobject$tpred<-object$predtime
predobject$alpha<-alpha
predobject
}

predict.ddpcomp<-function(object,newdata,alpha=0.05,tpred=NULL,...){
if(is.null(tpred)){
tpred<-object$predtime
}
tpred<-tpred/object$high.pct*10

xpred<-model.matrix(object$Terms,newdata,xlev = object$xlevels)
xpred<-xpred[,-1]

xpred<-(xpred-matrix(rep(object$xmean,times=nrow(xpred)),nrow=nrow(xpred), byrow=TRUE))/matrix(rep(2*object$xsd,times=nrow(xpred)),nrow=nrow(xpred), byrow=TRUE)	
predobject<-.Call('DPWeibull_predcompreg', PACKAGE = 'DPWeibull',
object$alpharec1,object$lambdascaled1,object$betarec1,
object$alpharec2,object$lambdascaled2,object$betarec2,object$prec,
xpred,tpred,alpha)
class(predobject)<-"predddpcomp"
predobject$tpred<-object$predtime
predobject$alpha<-alpha
predobject
}
