dpm<-function(tl,tr,event,high.pct,tpred,
                  burnin,iteration,
                  alpha00,
                  alpha0,
                  lambda00,
                  alphaalpha,alphalambda,
                  a,b,
                  addgroup,
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
  
  tl<-tl/high.pct*10
  tr<-tr/high.pct*10
  tpred<-tpred/high.pct*10
  nu<-rgamma(burnin+iteration+1,a,b)
  ngrp<-rep(1,burnin+iteration+1)
  result<-.Call('DPWeibull_noreg', PACKAGE = 'DPWeibull', 
        burnin, iteration, tl, tr, delta, pi,
        c, nm, alpha, lambda,
        lambda0, alpha00, alpha0, lambda00,
        alphaalpha, alphalambda,
        nu, ngrp, a, b, high.pct, tpred,addgroup, thin)
   result$Spred<-apply(result$S,2,median,na.rm=TRUE)
   result$Spredu<-apply(result$S,2,quantile,0.975,na.rm=TRUE)
   result$Spredl<-apply(result$S,2,quantile,0.025,na.rm=TRUE)
   result$dpred<-apply(result$d,2,median,na.rm=TRUE)
   result$dpredu<-apply(result$d,2,quantile,0.975,na.rm=TRUE)
   result$dpredl<-apply(result$d,2,quantile,0.025,na.rm=TRUE)
   result$hpred<-apply(result$h,2,median,na.rm=TRUE)
   result$hpredu<-apply(result$h,2,quantile,0.975,na.rm=TRUE)
   result$hpredl<-apply(result$h,2,quantile,0.025,na.rm=TRUE)
   result$tl<-tl/10*high.pct
   result$tr<-tr/10*high.pct
   result$pi<-pi
   result$delta<-delta
   result
  
}

