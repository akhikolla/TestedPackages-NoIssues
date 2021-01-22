compreg<-function(t,event,x,high.pct,predtime,indicator,
                  burnin,iteration,
                  alpha00,
                  alpha0,
                  lambda00,
                  alphaalpha,alphalambda,
                  a,b,
                  gamma0,gamma1,
                  addgroup,
                  thin, betasl){
  if(is.vector(x)){
    x<-matrix(x,ncol=1)
  }
  npts<-length(t)
  c<-rep(1,npts)
  nm<-c(npts,rep(0,npts-1))
  lambda01<-rgamma(npts,alpha00,lambda00)
  lambda1<-1:npts
  alpha1<-1:npts
  for(i in 1:npts){
    lambda1[i]<-rgamma(1,alpha0,lambda01[i])
    base<-ifelse(lambda1[i]==0,80,log(-log(0.05)/lambda1[i],base=25))
    if(base<80){
      alpha1[i]<-rtrunc(1, spec="gamma", a = max(base,0), b =80, 
                        shape=alphaalpha,rate=alphalambda)
    }else{
      alpha1[i]<-80
    }
  }
  lambda02<-rgamma(npts,alpha00,lambda00)
  lambda2<-1:npts
  alpha2<-1:npts
  for(i in 1:npts){
    lambda2[i]<-rgamma(1,alpha0,lambda02[i])
    base<-ifelse(lambda2[i]==0,80,log(-log(0.05)/lambda2[i],base=25))
    if(base<80){
      alpha2[i]<-rtrunc(1, spec="gamma", a = max(base,0), b =80, 
                        shape=alphaalpha,rate=alphalambda)
    }else{
      alpha2[i]<-80
    }
  }
  p<-rbeta(npts,1,1)
  beta1<-rep(0,npts*ncol(x))
  for(i in 1:npts*ncol(x)){
    beta1[i]<-rtrunc(1, spec="cauchy", a = -10, b =10, 
                    location=0.0,scale=betasl)
  }
  beta1<-matrix(beta1,ncol=ncol(x))
  
  beta2<-rep(0,npts*ncol(x))
  for(i in 1:npts*ncol(x)){
    beta2[i]<-rtrunc(1, spec="cauchy", a = -10, b =10, 
                    location=0.0,scale=betasl)
  }
  beta2<-matrix(beta2,ncol=ncol(x))
  
  t<-t/high.pct*10
  predtime<-predtime/high.pct*10
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
  result<-.Call('DPWeibull_compreg', PACKAGE = 'DPWeibull', burnin,iteration,
             t, x,event,
             c,nm,
             alpha1, lambda1, lambda01,
             alpha2, lambda2, lambda02,
             p, beta1, beta2,
             alpha00, alpha0, lambda00,alphaalpha, alphalambda,
             gamma0,gamma1, betasl,
             nu,ngrp,
             a, b,
             high.pct,predtime,
             addgroup,thin,xpred1, xpred2)
  xscale<-matrix(rep(2*xsd,length(predtime)),nrow=ncol(x))
  result$loghr.est<-matrix(apply(result$loghr,2,median,na.rm=TRUE),byrow=TRUE,nrow=ncol(x))/xscale
  result$loghrl<-matrix(apply(result$loghr,2,quantile,0.025,na.rm=TRUE),byrow=TRUE,nrow=ncol(x))/xscale
  result$loghru<-matrix(apply(result$loghr,2,quantile,0.975,na.rm=TRUE),byrow=TRUE,nrow=ncol(x))/xscale
  result$t<-t/10*high.pct
  result$xmean<-xmean
  result$xsd<-xsd
  result$xscale<-xscale
  result$event<-event
  result$indicator<-indicator
  result$covnames<-covnames
  result  
}

