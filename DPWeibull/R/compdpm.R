compdpm<-function(t,event,high.pct,times,
                  burnin,iteration,
                  alpha00,
                  alpha0,
                  lambda00,
                  alphaalpha,alphalambda,
                  a,b,
                  gamma0,gamma1,
                  addgroup,
                  thin){
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
  t<-t/high.pct*10
  times<-times/high.pct*10
  nu<-rgamma(burnin+iteration+1,a,b)
  ngrp<-rep(1,burnin+iteration+1)
 
  result<-.Call('DPWeibull_compnoreg', PACKAGE = 'DPWeibull', burnin,iteration,
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
             addgroup,thin)
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
   result$t<-t/10*high.pct
   result$event<-event
   result
}
