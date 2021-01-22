
reg<-function(tl,tr,event,x,high.pct,tpred,indicator,
                  burnin,iteration,
                  alpha00,
                  alpha0,
                  lambda00,
                  alphaalpha,alphalambda,
                  a,b,
                  addgroup,betasl,
                  thin){
  pi<-as.numeric(event==0)
  delta<-as.numeric((event==1)&(tl==tr))
  tr<-ifelse(event==3,tr,tl)
  npts<-length(tl)
  c<-rep(1,npts)
  nm<-c(npts,rep(0,npts-1))
  lambda0<-rgamma(npts,alpha00,lambda00)
  lambda<-1:npts
  alpha<-1:npts
  for(i in 1:npts){
    lambda[i]<-rgamma(1,alpha0,lambda0[i])
    base<-ifelse(lambda[i]==0,80,log(-log(0.05)/lambda[i],base=25))
    if(base<80){
      alpha[i]<-rtrunc(1, spec="gamma", a = max(base,0), b =80, 
                       shape=alphaalpha,rate=alphalambda)
    }else{
      alpha[i]<-80
    }
  }
  if(is.vector(x)){
    x<-matrix(x,ncol=1)
  }
  beta<-rtrunc(npts*ncol(x), spec="cauchy", a = -10, b =10, location=0.0,scale=betasl)
  beta<-matrix(beta,ncol=ncol(x))
  
  tl<-tl/high.pct*10
  tr<-tr/high.pct*10
  tpred<-tpred/high.pct*10
  xmean<-apply(x,2,mean)
  xsd<-ifelse(indicator,rep(0.5,length(xmean)),apply(x,2,sd))
  xpred1<-rep(1,length(xmean))
  xpred2<-rep(0,length(xmean))
  xpred1<-(xpred1-xmean)/2/xsd
  xpred2<-(xpred2-xmean)/2/xsd
  covnames<-names(xsd)
  x<-(x-matrix(rep(xmean,times=nrow(x)),nrow=nrow(x), byrow=TRUE))/matrix(rep(2*xsd,times=nrow(x)),nrow=nrow(x), byrow=TRUE)
  nu<-rgamma(burnin+iteration+1,a,b)
  ngrp<-rep(1,burnin+iteration+1)
  result<-.Call('DPWeibull_reg', PACKAGE = 'DPWeibull', 
        burnin, iteration, tl, tr, delta, pi,
        x, c, nm, alpha, lambda, beta,
        lambda0, alpha00, alpha0, lambda00,
        alphaalpha, alphalambda,
        nu, ngrp, a, b, high.pct, addgroup, betasl, xpred1,xpred2, tpred, thin)
  xscale<-matrix(rep(2*xsd,length(tpred)),nrow=ncol(x))
  result$loghr.est<-matrix(apply(result$loghr,2,median,na.rm=TRUE),byrow=TRUE,nrow=ncol(x))/xscale
  result$loghrl<-matrix(apply(result$loghr,2,quantile,0.025,na.rm=TRUE),byrow=TRUE,nrow=ncol(x))/xscale
  result$loghru<-matrix(apply(result$loghr,2,quantile,0.975,na.rm=TRUE),byrow=TRUE,nrow=ncol(x))/xscale
   result$tl<-tl/10*high.pct
   result$tr<-tr/10*high.pct
   result$pi<-pi
   result$delta<-delta
  result$xmean<-xmean
  result$xsd<-xsd
  result$xscale<-xscale
  result$indicator<-indicator
  result$covnames<-covnames
  result  
}
