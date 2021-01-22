onlymeanN = function(lower=rep(-Inf,length(mu)),upper=rep(Inf,length(mu)),mu,Sigma){
  p = length(mu)
  if(p==1){
    out = onlymeanNuni(a = lower,b = upper,mu = mu,Sigma = Sigma)  #OK
    return(out)
  }
  if(all(is.infinite(c(lower,upper)))){
    return(list(mean = mu))
  }
  bool = is.infinite(lower) & is.infinite(upper)
  #if exists (-Inf,Inf) limits
  if(sum(bool)>0){
    out = withinfs_mean(lower,upper,mu,Sigma,bool)
  }else{
    out = Vaida.LRIC.onlymean(a = lower,b = upper,mu = mu,Sigma = Sigma)
  }
  return(out)
}



###############################################################################################
###############################################################################################
#This function is optimized for the case when it DOES exist infinite values in the lower/upper
#truncation limits
###############################################################################################
###############################################################################################

Vaida.LRIC.onlymean<-function(a=rep(-Inf,length(mu)),b=rep(Inf,length(mu)),mu,Sigma){
  p=length(mu)
  if(p==1){
    return(onlymeanNuni(a,b,mu,Sigma))
  }else{
    #####
    Sigma = sym.matrix(Sigma)
    #####
    a1 <- (a-mu)/sqrt(diag(Sigma))
    b1 <- (b-mu)/sqrt(diag(Sigma))
    R <-  diag(1/sqrt(diag(Sigma)))%*%Sigma%*%diag(1/sqrt(diag(Sigma)))
    
    logp <- pmvnormt(lower = as.vector(a1),upper=as.vector(b1),sigma=R,uselog2 = TRUE)
    prob = 2^logp
    
    if(prob < 1e-100){
      #print("LRIC.Vaida corrector applied \n")
      return(corrector_onlymean(lower = a,upper = b,mu = mu,Sigma = Sigma,bw=36))
    }
    #mean
    qq = qfun(a1,b1,R)
    da = as.vector(qq$qa)
    db = as.vector(qq$qb)
    #EX <- -R%*%(da - db)/2^logp   # a vector with a length of p
    EX <- -R%*%log2ratio(da - db,logp)
    Eycens <- -diag(sqrt(diag(Sigma)))%*%EX+mu
    if(max(abs(Eycens))> 10*max(abs(c(a,b)[is.finite(c(a,b))])) | any(Eycens < a | Eycens > b)){
      #print("LRIC.Vaida corrector applied \n")
      return(corrector_onlymean(lower = a,upper = b,mu = mu,Sigma = Sigma,bw=36))
    }
  }
  return(list(mean=Eycens))
}
