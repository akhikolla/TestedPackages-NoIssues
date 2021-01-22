continue<-function(previous,iteration=1000,...) UseMethod("continue")

continue.dpm<-function(previous,...){
    continuedpm(previous,...)    
}

continuedpm<-function(previous,alpha= 0.05, simultaneous=FALSE,
		burnin=0,iteration=1000,
                  alpha00=1.354028,
                  alpha0=0.03501257,
                  lambda00=7.181247,
                  alphaalpha=0.2,alphalambda=0.1,
                  a=1,b=1,
                  addgroup=2,
                  thin=10){
	emptybasket<-previous$emptybasket
	allbaskets<-previous$allbaskets
	tl<-previous$tl
	tr<-previous$tr
	pi<-previous$pi
	delta<-previous$delta
	high.pct<-previous$high.pct
	tpred<-previous$predtime
	c<-previous$c
	nm<-previous$nm
	lastrow<-nrow(previous$alpharec)
	internalalpha<-alpha	
	alpha<-previous$alpharec[lastrow,]
	lambda<-previous$lambdascaled[lastrow,]
	lambda0<-previous$lambda0rec[lastrow,]
	ngrp<-previous$ngrp[length(previous$ngrp)]
	  npts<-length(tl)
	  tl<-tl/high.pct*10
	  tr<-tr/high.pct*10
	  tpred<-tpred/high.pct*10
	  nu<-rgamma(burnin+iteration+1,a,b)
	  ngrp<-c(ngrp, rep(1,burnin+iteration))
	  result<-.Call('DPWeibull_noreg_resume', PACKAGE = 'DPWeibull', 
		burnin, iteration, tl, tr, delta, pi,
		c, nm, alpha, lambda,
		lambda0, alpha00, alpha0, lambda00,
		alphaalpha, alphalambda,
		nu, ngrp, a, b, high.pct, tpred,addgroup, thin, emptybasket, allbaskets)
	   result$Spred<-apply(result$S,2,median,na.rm=TRUE)
	   result$Spredu<-apply(result$S,2,quantile,0.975,na.rm=TRUE)
	   result$Spredl<-apply(result$S,2,quantile,0.025,na.rm=TRUE)
	   result$dpred<-apply(result$d,2,median,na.rm=TRUE)
	   result$dpredu<-apply(result$d,2,quantile,0.975,na.rm=TRUE)
	   result$dpredl<-apply(result$d,2,quantile,0.025,na.rm=TRUE)
	   result$hpred<-apply(result$h,2,median,na.rm=TRUE)
	   result$hpredu<-apply(result$h,2,quantile,0.975,na.rm=TRUE)
	   result$hpredl<-apply(result$h,2,quantile,0.025,na.rm=TRUE)
	   result$predtime<-previous$predtime
	   result$tl<-previous$tl
	   result$tr<-previous$tr
	   result$pi<-pi
	   result$delta<-delta
	   result$high.pct<-high.pct
	   class(result)<-"dpm"
	   result$usertime<-previous$usertime
	   result$alpha<-internalalpha
	   result$simultaneous<-simultaneous
	   if(internalalpha!=0.05){
		result<-dpmdiffalpha(internalalpha,result)
		}
		if(simultaneous==TRUE){
		result$Sbandl<-confband(internalalpha,result$S)[1,]
		result$Sbandu<-confband(internalalpha,result$S)[2,]
		result$dbandl<-confband(internalalpha,result$d)[1,]
		result$dbandu<-confband(internalalpha,result$d)[2,]
		result$hbandl<-confband(internalalpha,result$h)[1,]
		result$hbandu<-confband(internalalpha,result$h)[2,]
		}  

   result
}
continue.ddp<-function(previous,...){
    continueddp(previous,...)    
}

continueddp<-function(previous,
	alpha= 0.05, simultaneous=FALSE,
	burnin=0,iteration=1000,
                  alpha00=1.354028,
                  alpha0=0.03501257,
                  lambda00=7.181247,
                  alphaalpha=0.2,alphalambda=0.1,
                  a=1,b=1,
                  addgroup=2,betasl=2.5,
                  thin=10){
	emptybasket<-previous$emptybasket
	allbaskets<-previous$allbaskets
	tl<-previous$tl
	tr<-previous$tr
	pi<-previous$pi
	delta<-previous$delta
	high.pct<-previous$high.pct
	tpred<-previous$predtime
	indicator<-previous$indicator
	xmean<-previous$xmean
	xsd<-previous$xsd
	c<-previous$c
	nm<-previous$nm
	lastrow<-nrow(previous$alpharec)
	internalalpha<-alpha	
	alpha<-previous$alpharec[lastrow,]
	lambda<-previous$lambdascaled[lastrow,]
	lambda0<-previous$lambda0rec[lastrow,]
	beta<-matrix(previous$betarec[lastrow,],byrow=TRUE,ncol=length(xmean))
	npts<-length(tl)
	tl<-tl/high.pct*10
	tr<-tr/high.pct*10
	tpred<-tpred/high.pct*10
	xpred1<-rep(1,length(xmean))
	xpred2<-rep(0,length(xmean))
	xpred1<-(xpred1-xmean)/2/xsd
	xpred2<-(xpred2-xmean)/2/xsd
  	x<-(previous$x-matrix(rep(xmean,times=nrow(previous$x)),nrow=nrow(previous$x), byrow=TRUE))/matrix(rep(2*xsd,times=nrow(previous$x)),nrow=nrow(previous$x), byrow=TRUE)
	nu<-rgamma(burnin+iteration+1,a,b)
	ngrp<-c(previous$ngrp, rep(1,burnin+iteration))
	result<-.Call('DPWeibull_reg_resume', PACKAGE = 'DPWeibull', 
	burnin, iteration, tl, tr, delta, pi,
	x, c, nm, alpha, lambda, beta,
	lambda0, alpha00, alpha0, lambda00,
	alphaalpha, alphalambda,
	nu, ngrp, a, b, high.pct, addgroup, betasl,
	xpred1,xpred2, tpred, thin,emptybasket, allbaskets)
	xscale<-matrix(rep(2*xsd,length(tpred)),nrow=ncol(x))			
	result$loghr.est<-matrix(apply(result$loghr,2,median,na.rm=TRUE),byrow=TRUE,nrow=ncol(x))/xscale
	result$loghrl<-matrix(apply(result$loghr,2,quantile,0.025,na.rm=TRUE),byrow=TRUE,nrow=ncol(x))/xscale
	result$loghru<-matrix(apply(result$loghr,2,quantile,0.975,na.rm=TRUE),byrow=TRUE,nrow=ncol(x))/xscale
	result$predtime<-previous$predtime
	result$tl<-previous$tl
	result$tr<-previous$tr
	result$pi<-pi
	result$delta<-delta
	result$high.pct<-high.pct
	result$xmean<-xmean
	result$xsd<-xsd
	result$x<-previous$x
	result$xscale<-previous$xscale
	result$covnames<-previous$covnames
	result$indicator<-indicator
	class(result)<-"ddp"
	result$usertime<-previous$usertime
	result$alpha<-internalalpha
	result$simultaneous<-simultaneous
	if(internalalpha!=0.05){
		result<-ddpdiffalpha(internalalpha,result)
	}
	if(simultaneous==TRUE){
		result$loghrbandl<-matrix(confband(internalalpha,result$loghr)[1,],byrow=TRUE,nrow=ncol(result$x))/result$xscale
		result$loghrbandu<-matrix(confband(internalalpha,result$loghr)[2,],byrow=TRUE,nrow=ncol(result$x))/result$xscale
	}
	result  
}

continue.dpmcomp<-function(previous,...){
    continuedpmcomp(previous,...)    
}

continuedpmcomp<-function(previous,alpha= 0.05, simultaneous=FALSE,
                  burnin=0,iteration=1000,
                  alpha00=1.354028,
                  alpha0=0.03501257,
                  lambda00=7.181247,
                  alphaalpha=0.2,alphalambda=0.1,
                  a=1,b=1,
                  gamma0=1,gamma1=1,
                  addgroup=2,
                  thin=1){
	emptybasket<-previous$emptybasket
	allbaskets<-previous$allbaskets
	t<-previous$t
	event<-previous$event
	high.pct<-previous$high.pct
	times<-previous$predtime
	c<-previous$c
	nm<-previous$nm
	lastrow<-nrow(previous$alpharec1)
	alpha1<-previous$alpharec1[lastrow,]
	lambda1<-previous$lambdascaled1[lastrow,]
	lambda01<-previous$lambda0rec1[lastrow,]
	alpha2<-previous$alpharec2[lastrow,]
	lambda2<-previous$lambdascaled2[lastrow,]
	lambda02<-previous$lambda0rec2[lastrow,]
	p<-previous$prec[lastrow,]
	ngrp<-previous$ngrp[length(previous$ngrp)]
        npts<-length(t)
  	t<-t/high.pct*10
  	times<-times/high.pct*10
  	nu<-rgamma(burnin+iteration+1,a,b)
  	ngrp<-c(ngrp, rep(1,burnin+iteration))
  	result<-.Call('DPWeibull_compnoreg_resume', PACKAGE = 'DPWeibull', burnin,iteration,
             t,event,
             c,nm,
             alpha1, lambda1, lambda01,
             alpha2, lambda2, lambda02,
             p,
             alpha00, alpha0, lambda00,alphaalpha, alphalambda,
             gamma0,gamma1,
             nu,ngrp,
             a, b,
             high.pct,times,
             addgroup,thin,emptybasket, allbaskets)
	   result$CIF1.est<-apply(result$CIF1,2,median,na.rm=TRUE)
	   result$CIF1u<-apply(result$CIF1,2,quantile,0.975,na.rm=TRUE)
	   result$CIF1l<-apply(result$CIF1,2,quantile,0.025,na.rm=TRUE)
	   result$CIF2.est<-apply(result$CIF2,2,median,na.rm=TRUE)
	   result$CIF2u<-apply(result$CIF2,2,quantile,0.975,na.rm=TRUE)
	   result$CIF2l<-apply(result$CIF2,2,quantile,0.025,na.rm=TRUE)
	   result$d1.est<-apply(result$d1,2,median,na.rm=TRUE)
	   result$d1u<-apply(result$d1,2,quantile,0.975,na.rm=TRUE)
	   result$d1l<-apply(result$d1,2,quantile,0.025,na.rm=TRUE)
	   result$d2.est<-apply(result$d2,2,median,na.rm=TRUE)
	   result$d2u<-apply(result$d2,2,quantile,0.975,na.rm=TRUE)
	   result$d2l<-apply(result$d2,2,quantile,0.025,na.rm=TRUE)
	   result$h1.est<-apply(result$h1,2,median,na.rm=TRUE)
	   result$h1u<-apply(result$h1,2,quantile,0.975,na.rm=TRUE)
	   result$h1l<-apply(result$h1,2,quantile,0.025,na.rm=TRUE)
	   result$h2.est<-apply(result$h2,2,median,na.rm=TRUE)
	   result$h2u<-apply(result$h2,2,quantile,0.975,na.rm=TRUE)
	   result$h2l<-apply(result$h2,2,quantile,0.025,na.rm=TRUE)
	   result$predtime<-previous$predtime
	   result$t<-t/10*high.pct
	   result$event<-event
	   result$high.pct<-high.pct
	   class(result)<-"dpmcomp"
	   result$usertime<-previous$usertime
	   result$alpha<-alpha
   	   result$simultaneous<-simultaneous        
	if(alpha!=0.05){
	result<-dpmcompdiffalpha(alpha,result)
	}
        if(simultaneous==TRUE){
        result$CIF1bandl<-confband(alpha,result$CIF1)[1,]
        result$CIF1bandu<-confband(alpha,result$CIF1)[2,]
        result$d1bandl<-confband(alpha,result$d1)[1,]
        result$d1bandu<-confband(alpha,result$d1)[2,]
        result$h1bandl<-confband(alpha,result$h1)[1,]
        result$h1bandu<-confband(alpha,result$h1)[2,]
        result$CIF2bandl<-confband(alpha,result$CIF2)[1,]
        result$CIF2bandu<-confband(alpha,result$CIF2)[2,]
        result$d2bandl<-confband(alpha,result$d2)[1,]
        result$d2bandu<-confband(alpha,result$d2)[2,]
        result$h2bandl<-confband(alpha,result$h2)[1,]
        result$h2bandu<-confband(alpha,result$h2)[2,]
	}  
   result
}

continue.ddpcomp<-function(previous,...){
    continueddpcomp(previous,...)    
}

continueddpcomp<-function(previous,alpha= 0.05, simultaneous=FALSE,
                  burnin=0,iteration=1000,
                  alpha00=1.354028,
                  alpha0=0.03501257,
                  lambda00=7.181247,
                  alphaalpha=0.2,alphalambda=0.1,
                  a=1,b=1,
                  gamma0=1,gamma1=1,
                  addgroup=2,
                  thin=10, betasl=2.5){
	emptybasket<-previous$emptybasket
	allbaskets<-previous$allbaskets
	t<-previous$t
	event<-previous$event
	high.pct<-previous$high.pct
	times<-previous$predtime
	indicator<-previous$indicator
	xmean<-previous$xmean
	xsd<-previous$xsd
	c<-previous$c
	nm<-previous$nm
       x<-(previous$x-matrix(rep(xmean,times=nrow(previous$x)),nrow=nrow(previous$x), byrow=TRUE))/matrix(rep(2*xsd,times=nrow(previous$x)),nrow=nrow(previous$x), byrow=TRUE)
	lastrow<-nrow(previous$alpharec1)
	alpha1<-previous$alpharec1[lastrow,]
	lambda1<-previous$lambdascaled1[lastrow,]
	lambda01<-previous$lambda0rec1[lastrow,]
	beta1<-matrix(previous$betarec1[lastrow,],byrow=TRUE,ncol=length(xmean))
	alpha2<-previous$alpharec2[lastrow,]
	lambda2<-previous$lambdascaled2[lastrow,]
	lambda02<-previous$lambda0rec2[lastrow,]
	beta2<-matrix(previous$betarec2[lastrow,],byrow=TRUE,ncol=length(xmean))
	p<-previous$prec[lastrow,]
	ngrp<-previous$ngrp[length(previous$ngrp)]
	  npts<-length(t)
	  t<-t/high.pct*10
	  times<-times/high.pct*10
	  xpred1<-rep(1,length(xmean))
	  xpred2<-rep(0,length(xmean))
	  xpred1<-(xpred1-xmean)/2/xsd
	  xpred2<-(xpred2-xmean)/2/xsd
	  nu<-rgamma(burnin+iteration+1,a,b)
	  ngrp<-c(ngrp, rep(1,burnin+iteration))
 
	  result<-.Call('DPWeibull_compreg_resume', PACKAGE = 'DPWeibull', burnin,iteration,
             t, x,event,
             c,nm,
             alpha1, lambda1, lambda01,
             alpha2, lambda2, lambda02,
             p, beta1, beta2,
             alpha00, alpha0, lambda00,alphaalpha, alphalambda,
             gamma0,gamma1, betasl,
             nu,ngrp,
             a, b,
             high.pct,times,
             addgroup,thin,xpred1, xpred2,emptybasket, allbaskets)
  xscale<-matrix(rep(2*xsd,length(times)),nrow=ncol(x))
  result$loghr.est<-matrix(apply(result$loghr,2,median,na.rm=TRUE),byrow=TRUE,nrow=ncol(x))/xscale
  result$loghrl<-matrix(apply(result$loghr,2,quantile,0.025,na.rm=TRUE),byrow=TRUE,nrow=ncol(x))/xscale
  result$loghru<-matrix(apply(result$loghr,2,quantile,0.975,na.rm=TRUE),byrow=TRUE,nrow=ncol(x))/xscale
  result$predtime<-previous$predtime
  result$t<-t/10*high.pct
  result$high.pct<-high.pct
  result$xmean<-xmean
  result$xsd<-xsd
  result$x<-previous$x
  result$xscale<-previous$xscale
  result$covnames<-previous$covnames
   result$event<-event
  result$indicator<-indicator
   class(result)<-"ddpcomp"
   result$usertime<-previous$usertime

      result$alpha<-alpha
   result$simultaneous<-simultaneous 
        if(alpha!=0.05){
	result<-ddpdiffalpha(alpha,result)
	}
        if(simultaneous==TRUE){
        result$loghrbandl<-matrix(confband(alpha,result$loghr)[1,],byrow=TRUE,nrow=ncol(result$x))/result$xscale
  	result$loghrbandu<-matrix(confband(alpha,result$loghr)[2,],byrow=TRUE,nrow=ncol(result$x))/result$xscale
	}
  result  
}

