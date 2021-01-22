onlymeanESN = function(lower = rep(-Inf,length(mu)),upper = rep(Inf,length(mu)),mu = mu,Sigma = Sigma,lambda = lambda,tau = tau){
  p = length(mu)
  if(p==1){
    out = onlymeanESNuni(lower,upper,mu,Sigma,lambda,tau)
    return(out)
  }
  #tau goes to infinite
  tautil<-tau/sqrt(1+sum(lambda^2))
  if(tautil< -36){
    #print("normal aproximation")
    Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
    return(onlymeanN(lower = lower,upper = upper,mu = mu - tautil*Delta,Sigma = Sigma - Delta%*%t(Delta)))
  }
  if(all(is.infinite(lower))){
    if(all(is.infinite(upper))){
      #No truncating at all
      return(ESN.NOTRUNC.onlymean(mu,Sigma,lambda,tau))  #OK
    }else
    {
      #Right censoring
      SS   = sqrtm(Sigma)
      tautil = tau/sqrt(1+sum(lambda^2))
      #xi = pnorm(tautil)
      varpsi = lambda/sqrt(1+sum(lambda^2))
      Omega  = cbind(rbind(Sigma,-t(varpsi)%*%SS),rbind(-SS%*%varpsi,1))
      rownames(Omega) <- colnames(Omega)

      bool = is.infinite(upper)
      #if exists (-Inf,Inf) limits
      if(sum(bool)>0){
        out = withinfs_mean(upper = c(upper,tautil),mu = c(mu,0),Sigma = Omega,bool = c(bool,FALSE))
      }else
      {
        out = Vaida.LRIC.onlymean(b = c(upper,tautil),mu = c(mu,0),Sigma = Omega)
      }
    }
  }else
  {
    if(all(is.infinite(upper))){
      #Left censoring
      SS   = sqrtm(Sigma)
      tautil = tau/sqrt(1+sum(lambda^2))
      #xi = pnorm(tautil)
      varpsi = lambda/sqrt(1+sum(lambda^2))
      Omega  = cbind(rbind(Sigma,t(varpsi)%*%SS),rbind(SS%*%varpsi,1))
      rownames(Omega) <- colnames(Omega)

      bool = is.infinite(lower)
      #if exists (-Inf,Inf) limits
      if(sum(bool)>0){
        out = withinfs_mean(upper = c(-lower,tautil),mu = c(-mu,0),Sigma = Omega,bool = c(bool,FALSE))
      }else
      {
        out = Vaida.LRIC.onlymean(b = c(-lower,tautil),mu = c(-mu,0),Sigma = Omega)
      }
      out$mean = -out$mean
    }else
    {
      SS   = sqrtm(Sigma)
      tautil = tau/sqrt(1+sum(lambda^2))
      #xi = pnorm(tautil)
      varpsi = lambda/sqrt(1+sum(lambda^2))
      Omega  = cbind(rbind(Sigma,-t(varpsi)%*%SS),rbind(-SS%*%varpsi,1))
      rownames(Omega) <- colnames(Omega)
      #intervalar censoring
      if(all(is.finite(c(lower,upper)))){
        #no infinites #all intervalar truncated
        out = Vaida.LRIC.onlymean(a = c(lower,-10^7),b = c(upper,tautil),mu = c(mu,0),Sigma = Omega)
      }else
      {
        #All kind of censoring
        bool = is.infinite(lower) & is.infinite(upper)
        #if exists (-Inf,Inf) limits
        if(sum(bool)>0){
          out = withinfs_mean(c(lower,-Inf),c(upper,tautil),c(mu,0),Omega,bool = c(bool,FALSE))
        }else{
          out = Vaida.LRIC.onlymean(a = c(lower,-Inf),b = c(upper,tautil),mu = c(mu,0),Sigma = Omega)
        }
      }
    }
  }
  return(list(mean = matrix(out$mean[-(p+1)],p,1)))
}
