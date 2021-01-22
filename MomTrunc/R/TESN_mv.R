###############################################################################################
###############################################################################################
#This function compute the mean and the var-cov matrix for a TESN distribution using the
#Normal augmented model
###############################################################################################
###############################################################################################

meanvarESN7 = function(lower=rep(-Inf,length(mu)),upper=rep(Inf,length(mu)),mu,Sigma,lambda,tau){
  p = length(mu)
  if(p==1){
    out = meanvarESNuni(a = lower,b = upper,mu = mu,Sigma = Sigma,lambda = lambda,tau = tau)  #OK
    return(out)
  }
  #tau goes to infinite
  tautil<-tau/sqrt(1+sum(lambda^2))
  if(tautil< -36){
    #print("normal aproximation")
    Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
    return(meanvarN7(lower = lower,upper = upper,mu = mu - tautil*Delta,Sigma = Sigma - Delta%*%t(Delta)))
  }
  if(all(is.infinite(lower))){
    if(all(is.infinite(upper))){
      #No truncating at all
      return(ESN.NOTRUNC(mu,Sigma,lambda,tau))  #OK
    }else{
      #Right censoring
      SS   = sqrtm(Sigma)
      varpsi = lambda/sqrt(1+sum(lambda^2))
      Omega  = cbind(rbind(Sigma,-t(varpsi)%*%SS),rbind(-SS%*%varpsi,1))
      rownames(Omega) <- colnames(Omega)


      bool = is.infinite(upper)
      #if exists (-Inf,Inf) limits
      if(sum(bool)>0){
        out = withinfs(upper = c(upper,tautil),mu = c(mu,0),Sigma = Omega,bool = c(bool,FALSE))
      }else{
        if(p<4){
          out = Kan.RC(b = c(upper,tautil),mu = c(mu,0),Sigma = Omega)
        }else{
          out = Vaida.RC(b = c(upper,tautil),mu = c(mu,0),Sigma = Omega)   #OK
        }
      }
    }
  }else{
    if(all(is.infinite(upper))){
      #Left censoring
      SS   = sqrtm(Sigma)
      varpsi = lambda/sqrt(1+sum(lambda^2))
      Omega  = cbind(rbind(Sigma,t(varpsi)%*%SS),rbind(SS%*%varpsi,1))
      rownames(Omega) <- colnames(Omega)

      bool = is.infinite(lower)
      #if exists (-Inf,Inf) limits
      if(sum(bool)>0){
        out = withinfs(upper = c(-lower,tautil),mu = c(-mu,0),Sigma = Omega,bool = c(bool,FALSE))
      }else{
        if(p<4){
          out = Kan.RC(b = c(-lower,tautil),mu = c(-mu,0),Sigma = Omega) #OK
        }else{
          out = Vaida.RC(b = c(-lower,tautil),mu = c(-mu,0),Sigma = Omega) #OK
        }
      }
      out$mean = -out$mean
    }else{
      SS   = sqrtm(Sigma)
      #tautil = tau/sqrt(1+sum(lambda^2))
      #xi = pnorm(tautil)
      varpsi = lambda/sqrt(1+sum(lambda^2))
      Omega  = cbind(rbind(Sigma,-t(varpsi)%*%SS),rbind(-SS%*%varpsi,1))
      rownames(Omega) <- colnames(Omega)
      #intervalar censoring
      if(all(is.finite(c(lower,upper)))){
        #no infinites #all intervalar truncated
        if(p<4){
          out = Kan.LRIC(a = c(lower,-Inf),b = c(upper,tautil),mu = c(mu,0),Sigma = Omega)
        }else{
          out = Vaida.LRIC(a = c(lower,-Inf),b = c(upper,tautil),mu = c(mu,0),Sigma = Omega)
        }
      }else
      {
        #All kind of censoring
        bool = is.infinite(lower) & is.infinite(upper)
        #if exists (-Inf,Inf) limits
        if(sum(bool)>0){
          out = withinfs(c(lower,-Inf),c(upper,tautil),c(mu,0),Omega,bool = c(bool,FALSE))
        }else{
          if(p<4){
            out = Kan.LRIC(a = c(lower,-Inf),b = c(upper,tautil),mu = c(mu,0),Sigma = Omega)
          }else{
            out = Vaida.LRIC(a = c(lower,-Inf),b = c(upper,tautil),mu = c(mu,0),Sigma = Omega)
          }
        }
      }
    }
  }
  return(list(mean = matrix(out$mean[-(p+1)],p,1),EYY = out$EYY[-(p+1),-(p+1)],varcov = out$varcov[-(p+1),-(p+1)]))
}
