onlymeanTall = function(lower=rep(-Inf,length(mu)),upper=rep(Inf,length(mu)),mu,Sigma,nu){
  p = length(mu)
  if(nu >= 3){
    return(onlymeanT(a = lower,b = upper,mu = mu,Sigma = Sigma,nu=nu))  #OK
  }else{
    if(all(is.infinite(lower))){
      if(all(is.infinite(upper))){
        #No truncating at all
        return(list(mean = mu))
      }else{
        #Right censoring
        return(utmvtmu(b = upper,mu = mu,S = Sigma,nu=nu))
      }
    }else{
      if(all(is.infinite(upper))){
        #Left censoring
        return(ltmvtmu(a = lower,mu = mu,S = Sigma,nu=nu))
      }
      else{
        return(dtmvtmu(a = lower,b = upper,mu = mu,S = Sigma,nu=nu))
      }
    }
  }
}
