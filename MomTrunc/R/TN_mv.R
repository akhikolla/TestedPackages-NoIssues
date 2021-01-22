###############################################################################################
###############################################################################################
#This function compute the mean and the var-cov matrix for a TN distribution
###############################################################################################
###############################################################################################

meanvarN7 = function(lower=rep(-Inf,length(mu)),upper=rep(Inf,length(mu)),mu,Sigma){
  p = length(mu)
  if(p==1){
    out = meanvarNuni(a = lower,b = upper,mu = mu,Sigma = Sigma)  #OK
    return(out)
  }
  if(p<10){
    #meanvarN
    if(all(is.infinite(lower))){
      if(all(is.infinite(upper))){
        #No truncating at all
        return(list(mean = mu,EYY = Sigma + mu%*%t(mu),varcov = Sigma))
      }else{
        #Right censoring
        bool = is.infinite(upper)
        #if exists (-Inf,Inf) limits
        if(sum(bool)>0){
          out = withinfs(upper = upper,mu = mu,Sigma = Sigma,bool = bool)
        }else{
          out = Kan.RC(b = upper,mu = mu,Sigma = Sigma)
        }
      }
    }else{
      if(all(is.infinite(upper))){
        #Left censoring
        bool = is.infinite(lower)
        #if exists (-Inf,Inf) limits
        if(sum(bool)>0){
          out = withinfs(upper = -lower,mu = -mu,Sigma = Sigma,bool = bool)
          out$mean = -out$mean
        }else{
          out = Kan.RC(b = -lower,mu = -mu,Sigma = Sigma) #OK
          out$mean = -out$mean
        }
      }else{
        #intervalar censoring
        if(all(is.finite(c(lower,upper)))){
          #no infinites #all intervalar truncated
          #print("IC")
          out = Kan.IC(a = lower,b = upper,mu = mu,Sigma = Sigma)
        }else{
          #All kind of censoring
          bool = is.infinite(lower) & is.infinite(upper)
          #if exists (-Inf,Inf) limits
          if(sum(bool)>0){
            out = withinfs(lower,upper,mu,Sigma,bool)
          }else{
            out = Kan.LRIC(a = lower,b = upper,mu = mu,Sigma = Sigma)
          }
        }
      }
    }
  }else{
    #vaida
    if(all(is.infinite(lower))){
      if(all(is.infinite(upper))){
        #No truncating at all
        return(list(mean = mu,EYY = Sigma + mu%*%t(mu),varcov = Sigma))
      }else{
        #Right censoring
        bool = is.infinite(upper)
        #if exists (-Inf,Inf) limits
        if(sum(bool)>0){
          out = withinfs(upper = upper,mu = mu,Sigma = Sigma,bool = bool)
        }else{
          out = Vaida.RC(b = upper,mu = mu,Sigma = Sigma)   #OK
        }
      }
    }else{
      if(all(is.infinite(upper))){
        #Left censoring
        bool = is.infinite(lower)
        #if exists (-Inf,Inf) limits
        if(sum(bool)>0){
          out = withinfs(upper = -lower,mu = -mu,Sigma = Sigma,bool = bool)
          out$mean = -out$mean
        }else{
          out = Vaida.RC(b = -lower,mu = -mu,Sigma = Sigma) #OK
          out$mean = -out$mean
        }
      }else{
        #intervalar censoring
        if(all(is.finite(c(lower,upper)))){
          #no infinites #all intervalar truncated
          out = Vaida.IC(a = lower,b = upper,mu = mu,Sigma = Sigma)
        }else{
          #All kind of censoring
          bool = is.infinite(lower) & is.infinite(upper)
          #if exists (-Inf,Inf) limits
          if(sum(bool)>0){
            out = withinfs(lower,upper,mu,Sigma,bool)
          }else{
            out = Vaida.LRIC(a = lower,b = upper,mu = mu,Sigma = Sigma)
          }
        }
      }
    }
  }
  return(out)
}
