dtmvtmu = function(a,b,mu,S,nu){

  n = length(mu)
  s = sqrt(diag(as.matrix(S)))
  a1 = a-mu
  b1 = b-mu
  as = a1/s
  bs = b1/s

  log2p = pmvt.genz(lower = a1,upper = b1,nu = nu,sigma = S,uselog2 = TRUE)$Estimation
  p = 2^log2p


  if(n == 1){
    muY  = 0
    if(is.infinite(a)||is.infinite(b)){
      if(nu<=1){
        muY = NaN
      }
    }
    kappa = log2ratio(s*gamma((nu+1)/2)/gamma(nu/2)/sqrt(nu*pi),log2p)


    if(is.infinite(a)){
      a  = 0
      av = 0
    }else{
      av = (1+as^2/nu)^(-(nu-1)/2)
    }
    if(is.infinite(b)){
      b =  0
      bv = 0
    }else{
      bv = (1+bs^2/nu)^(-(nu-1)/2)
    }
    if(!is.na(muY)){
      if(nu == 1){
        muY = mu+s*log((1+bs^2)/(1+as^2))/(atan(bs)-atan(as))/2
      }else{
        muY = mu+kappa*nu/(nu-1)*(av-bv)
      }
    }
    return(list(mean = muY))
  }

  ind = seq_len(n)[is.infinite(a)|is.infinite(b)]
  d = n-length(ind)+nu
  if(nu==d && nu<=1){
    return(list(mean = matrix(NaN,n,1)))
  }
  qout = qfunT(a1,b1,S,nu)
  qa = qout$qa
  qb = qout$qb
  q = qa-qb
  muY = mu+S%*%log2ratio(q,log2p)

  #Check existence of moments
  if(nu<=2){
    if(d<=1){
      muY[ind] = NaN
    }
  }

  return(list(mean = muY))
}

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

ltmvtmu = function(a,mu,S,nu){

  n = length(mu)

  if(nu<=1){
    return(list(mean = matrix(NaN,n,1)))
  }
  s  = sqrt(diag(as.matrix(S)))
  a1 = a-mu
  as = a1/s
  c  = sqrt(nu)*gamma((nu-1)/2)/(2*sqrt(pi)*gamma(nu/2))*(1+as^2/nu)^(-(nu-1)/2)
  if(n==1){
    if(is.infinite(a)){
      muY = mu
    }else{
      c1 = s*c/pt(-as,nu)
      muY = mu+c1
    }
    return(list(mean = muY))
  }

  log2p = pmvt.genz(lower = a1,nu = nu,sigma = S,uselog2 = TRUE)$Estimation
  #p = 2^log2p
  #p = pmvt(lower = a1,df = nu,sigma = S)[1]
  q = qfunT_a(a,mu,S,nu)
  #muY = mu+S%*%q/p
  muY = mu+S%*%log2ratio(q,log2p)

  return(list(mean = muY))
}


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

utmvtmu = function(b,mu,S,nu){
  out = ltmvtmuvar(-b,-mu,S,nu)
  out$muY = -out$muY
  return(out)
}
