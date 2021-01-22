dtmvtmuvar = function(a,b,mu,S,nu){
  
  S = as.matrix(S)
  n = length(mu)
  s = sqrt(diag(S))
  a1 = a-mu
  b1 = b-mu
  as = a1/s
  bs = b1/s
  
  log2p = pmvt.genz(lower = a1,upper = b1,nu = nu,sigma = as.matrix(S),uselog2 = TRUE,N = 799)$Estimation
  p = 2^log2p
  
  
  if(n == 1){
    muY  = 0
    varY = 0
    if(is.infinite(a)||is.infinite(b)){
      if(nu<=1){
        muY = NaN
        varY = NaN
      }else{
        if(nu <= 2){varY = NaN}
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
      if(!is.na(varY)){
        if (nu == 2){
          c = asinh(bs/sqrt(2))-asinh(as/sqrt(2))+(as*av-bs*bv)/sqrt(2)
          varY = log2ratio(c*S,log2p)-(muY-mu)^2
        }else{
          varY = ((nu-3)*mu*muY+(nu*S+mu^2)+nu*kappa*(a*av-b*bv))/(nu-2)-muY^2
        }
      }
    }
    return(list(mean = muY,EYY = varY +muY^2,varcov = varY))
  }
  
  ind = seq_len(n)[is.infinite(a)|is.infinite(b)]
  d = n-length(ind)+nu
  if(nu==d && nu<=1){
    return(list(mean = matrix(NaN,n,1),EYY = matrix(NaN,n,n),varcov = matrix(NaN,n,n)))
  }
  qout = qfunT(a1,b1,S,nu)
  qa = qout$qa
  qb = qout$qb
  q = qa-qb
  muY = mu+S%*%log2ratio(q,log2p)
  
  if(n == 2){
    rho = S[1,2]/(s[1]*s[2])
    if(nu==1 && d<=2){
      varY = matrix(NaN,2,2)
      # %  One pair of a and b is finite
      if(is.infinite(a[2]) && is.infinite(b[2])){  #Finite a(1) and b(1), a(2)=-Inf, b(2)=Inf
        varY[1,1] = mu[1]*(2*muY[1]-mu[1])+S[1,1]*((bs[1]-as[1])/pi/p-1)
        varY[1,2] = (mu[2]-S[1,2]*mu[1]/S[1,1])*muY[1]+S[1,2]/S[1,1]*varY[1,1]
      }else{
        if(is.infinite(a[1]) && is.infinite(b[1])){  #% Finite a(2) and b(2), a(1)=-Inf, b(1)=Inf
          varY[2,2] = mu[2]*(2*muY[2]-mu[2])+S[2,2]*((bs[2]-as[2])/pi/p-1)
          varY[1,2] = (mu[1]-S[1,2]*mu[2]/S[2,2])*muY[2]+S[1,2]/S[2,2]*varY[2,2]
        }else{
          if(is.infinite(a[2])||is.infinite(b[2])){  #% Finite a(1) and b(1)
            xl = as[1]
            xu = bs[1]
            if(is.infinite(a[2])){
              yy = bs[2]
            }else{
              yy = as[2]
            }
            J0 = (qa[2]-qb[2])*s[2]
          }else{                          # Finite a(2) and b(2)
            xl = as[2]
            xu = bs[2]
            if(is.infinite(a[1])){
              yy = bs[1]
            }else{
              yy = as[1]
            }
            J0 = (qa[1]-qb[1])*s[1]
          }
          ca = sqrt(1-rho^2+yy^2-2*rho*yy*xl+xl^2)
          cb = sqrt(1-rho^2+yy^2-2*rho*yy*xu+xu^2)
          J1 = (ca-cb)/(2*pi)
          if(is.infinite(a[1]) || is.infinite(a[2])){
            qq = (rho*J1-yy*J0*(1-rho^2)+(xu-xl)/(2*pi))/p-1
            qq1 = log2ratio(J1+yy*rho*J0,log2p)
          }else{
            qq = (-rho*J1-yy*J0*(1-rho^2)+(xu-xl)/(2*pi))/p-1
            qq1 = log2ratio(-J1+yy*rho*J0,log2p)
          }
          if(is.infinite(a[2])||is.infinite(b[2])){
            varY[1,1] = mu[1]*(2*muY[1]-mu[1])+S[1,1]*qq
          }else{
            varY[2,2] = mu[2]*(2*muY[2]-mu[2])+S[2,2]*qq
          }
          varY[1,2] = mu[1]*muY[2]+mu[2]*muY[1]-mu[1]*mu[2]+S[1,2]*qq+s[1]*s[2]*(1-rho^2)*qq1
        }
      }
      varY[2,1] = varY[1,2]
      #varY = varY-muY*t(muY)
      return(list(mean = muY, EYY = varY, varcov = varY-muY*t(muY)))
    }
    if(nu>2||(d>nu && nu != 2)){
      #  A faster way of computing nu/(nu-2)*mvtcdf(as/c,bs/c,R,nu-2), where c=sqrt(nu/(nu-2))
      ind1 = is.finite(a)
      ind2 = is.finite(b)
      p0 = (nu*p + sum(a1[ind1]*qa[ind1]) - sum(b1[ind2]*qb[ind2]))/(nu-2)
      
    }else{
      if(d>nu && nu==2){   #check if there exists a pair of (a_i,b_i) that is finite
        if(is.infinite(a[2])&&is.infinite(b[2])){
          p0 = asinh(bs[1]/sqrt(2))-asinh(as[1]/sqrt(2))
        }else{
          if(is.infinite(a[1]) && is.infinite(b[1])){
            p0 = asinh(bs[2]/sqrt(2))-asinh(as[2]/sqrt(2))
          }else{
            faux = function(w1){
              return(1/sqrt(2+w1^2)*
                       (atan((bs[2]-rho*w1)/sqrt(1-rho^2)/sqrt(2+w1^2)) -
                          atan((as[2]-rho*w1)/sqrt(1-rho^2)/sqrt(2+w1^2))))
            }
            p0 = integrate(faux,as[1],bs[1],abs.tol = 1e-10,rel.tol = 1e-10)$value/pi
          }
        }
      }else{
        p0 = 0
      }
    }
    Si = solve(S)
    c = c(t(a1)%*%Si%*%a1,
          t(c(a1[1],b1[2]))%*%Si%*%c(a1[1],b1[2]),
          t(c(b1[1],a1[2]))%*%Si%*%c(b1[1],a1[2]),
          t(b1)%*%Si%*%b1)
    if(nu==2){
      p1 = -log(2+c)
    }else{
      p1 = nu/(nu-2)*(1+c/nu)^(-nu/2+1)
    }
    p1[is.infinite(c)|is.na(c)] = 0
    h = sqrt(1-rho^2)*(p1[1]-p1[2]-p1[3]+p1[4])/(2*pi)
    W = matrix(0,2,2)
    a[is.infinite(a)]   = 0
    b[is.infinite(b)]   = 0
    a1[is.infinite(a1)] = 0
    b1[is.infinite(b1)] = 0
    #  W = mu*q'*S+S*(D-p0*eye(2))
    W[1,1] = S[1,1]*((a[1]+mu[1])*qa[1]-(b[1]+mu[1])*qb[1])+
      S[1,2]*((2*mu[1]+S[1,2]/S[2,2]*a1[2])*qa[2]-(2*mu[1]+S[1,2]/S[2,2]*b1[2])*qb[2])+S[1,1]*rho*h
    W[2,2] = S[2,2]*((a[2]+mu[2])*qa[2]-(b[2]+mu[2])*qb[2])+
      S[1,2]*((2*mu[2]+S[1,2]/S[1,1]*a1[1])*qa[1]-(2*mu[2]+S[1,2]/S[1,1]*b1[1])*qb[1])+S[2,2]*rho*h
    W[1,2] = S[1,1]*mu[2]*q[1]+S[2,2]*mu[1]*q[2]+
      S[1,2]*(a[1]*qa[1]+a[2]*qa[2]-b[1]*qb[1]-b[2]*qb[2])+s[1]*s[2]*h
    W[2,1] = W[1,2]
    varY = log2ratio(p0*S + W,log2p) + mu%*%t(mu) - muY%*%t(muY)
    
  }else{
    h = sqrt(nu)/(2*sqrt(pi)*gamma(nu/2))
    
    if(nu==1 && sum(is.infinite(a*b)) == n-1){
      # A special case that requires interpolation
      D = matrix(NaN,n,n)
      i = is.finite(a*b)
      nu1 = nu+1e-4
      nu2 = nu-1e-4
      # nu1 = nu + 1
      # nu2 = nu - 1
      qout1  = qfunT(a1,b1,S,nu1)#error
      qout2  = qfunT(a1,b1,S,nu2)#error
      # d1     = lfun(a1,b1,nu1*S,nu1-2)*nu1/(2*gamma(nu1/2))+a[i]*qout1$qa[i]-b[i]*qout1$qb[i]
      # d2     = lfun(a1,b1,nu2*S,nu2-2)*nu2/(2*gamma(nu2/2))+a[i]*qout2$qa[i]-b[i]*qout2$qb[i]
      d1     = lfun(a1,b1,nu1*S,nu1-2)*nu1/(2*gamma(nu1/2))+a[i]*qout1$qa[i]-b[i]*qout1$qb[i]
      d2     = lfun(a1,b1,nu2*S,nu2-2)*nu2/(2*gamma(nu2/2))+a[i]*qout2$qa[i]-b[i]*qout2$qb[i]
      D[i,i] = (d1+d2)/2
    }else{
      D = lfun(a1,b1,nu*S,nu-2)*nu/(2*gamma(nu/2))*diag(n)
      for(i in 1:n){
        if(is.finite(a[i])){
          D[i,i] = D[i,i]+a[i]*qa[i]
        }
        if(is.finite(b[i])){
          D[i,i] = D[i,i]-b[i]*qb[i]
        }
      }
    }
    
    for(i in 1:n){
      RR = S[-i,-i]-S[-i,i]%*%t(S[i,-i])/S[i,i]
      if(is.finite(a[i]*b[i]) && ((nu==2 && all(is.infinite(a[-i]*b[-i]))) || (nu==1 && sum(is.infinite(a[-i]*b[-i])) == n-2))){
        #Special cases that require interpolation
        nu1 = nu + 1e-1
        nu2 = nu - 1e-1
        # nu1 = nu + 1
        # nu2 = nu - 1
        ma = mu[-i]+S[-i,i]/S[i,i]*a1[i]
        mb = mu[-i]+S[-i,i]/S[i,i]*b1[i]
        wa1 = hfun(a[-i],b[-i],ma,RR*(nu1+as[i]^2),nu1-1)*(1+as[i]^2/nu1)^(-(nu1-1)/2)-
          hfun(a[-i],b[-i],mb,RR*(nu1+bs[i]^2),nu1-1)*(1+bs[i]^2/nu1)^(-(nu1-1)/2)
        wa2 = hfun(a[-i],b[-i],ma,RR*(nu2+as[i]^2),nu2-1)*(1+as[i]^2/nu2)^(-(nu2-1)/2)-
          hfun(a[-i],b[-i],mb,RR*(nu2+bs[i]^2),nu2-1)*(1+bs[i]^2/nu2)^(-(nu2-1)/2)
        wa = (wa1+wa2)*h/(2*s[i])
      }else{
        if(is.finite(a[i])){
          ma = mu[-i]+S[-i,i]/S[i,i]*a1[i]
          wa = hfun(a[-i],b[-i],ma,RR*(nu+as[i]^2),nu-1)*(1+as[i]^2/nu)^(-(nu-1)/2)*h/s[i]
        }else{
          wa = matrix(0,n-1,1)
        }
        if(is.finite(b[i])){
          mb = mu[-i]+S[-i,i]/S[i,i]*b1[i]
          wa = wa-hfun(a[-i],b[-i],mb,RR*(nu+bs[i]^2),nu-1)*(1+bs[i]^2/nu)^(-(nu-1)/2)*h/s[i]
        }
      }
      D[i,-i] = wa
    }
    
    varY = log2ratio(S%*%(D-q%*%t(muY)),log2p)
    #  Need to understand why this statement is needed
    if(nu == 1 && sum(is.infinite(a*b)) == n-1){
      i = is.finite(a*b)
      varY[i,] = varY[,i]
    }
  }
  
  varY = (varY + t(varY))/2
  
  #Check existence of moments
  if(nu<=2){
    if(d<=1){
      varY[ind,] = NaN
      varY[,ind] = NaN
      muY[ind]   = NaN
    }else{
      if(d<=2){
        varY[ind,ind] = NaN
      }
    }
  }
  
  return(list(mean = muY,EYY = varY + muY%*%t(muY),varcov = varY))
}

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

ltmvtmuvar = function(a,mu,S,nu){
  
  n = length(mu)
  
  if(nu<=1){
    return(list(mean = matrix(NaN,n,1),EYY = matrix(NaN,n,n),varcov = matrix(NaN,n,n)))
  }
  
  S = as.matrix(S)
  s = sqrt(diag(S))
  
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
    if(nu<=2){
      varY = NaN
    }else{
      if(is.infinite(a)){
        varY = nu/(nu-2)*S
      }else{
        varY = nu/(nu-2)*S+s*as*c1*(nu-1)/(nu-2)-c1^2
      }
    }
    return(list(mean = muY,EYY = varY +muY^2,varcov = varY))
  }
  
  log2p = pmvt.genz(lower = a1,nu = nu,sigma = S,uselog2 = TRUE,N = 799)$Estimation
  #p = 2^log2p
  #p = pmvt(lower = a1,df = nu,sigma = S)[1]
  q = qfunT_a(a,mu,S,nu)
  #muY = mu+S%*%q/p
  muY = mu+S%*%log2ratio(q,log2p)
  
  if(nu<=2){
    varY = matrix(NaN,n,n)
  }else{
    #  A faster way of computing nu/(nu-2)*mvtcdf(-sqrt((nu-2)/nu)*as,R,nu-2)
    ind2 = is.finite(a)
    p1 = as.numeric((nu*2^log2p+t(a[ind2]-mu[ind2])%*%q[ind2])/(nu-2))
    D = p1*diag(n)
    for(i in 1:n){
      if(is.finite(a[i])){
        D[i,i] = D[i,i]+a[i]*q[i]
        mu1 = mu[-i]+S[-i,i]*as[i]/s[i]
        S1 = (nu+as[i]^2)/(nu-1)*(S[-i,-i]-S[-i,i]%*%t(S[i,-i])/S[i,i])
        D[i,-i] = q[i]*mu1 + c[i]/s[i]*S1%*%qfunT_a(a = a[-i],mu = mu1,S = S1,nu-1)
      }
    }
    varY = log2ratio(S%*%(D-q%*%t(muY)),log2p)
  }
  
  varY = (varY + t(varY))/2
  
  return(list(mean = muY,EYY = varY + muY%*%t(muY),varcov = varY))
}


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

utmvtmuvar = function(b,mu,S,nu){
  out = ltmvtmuvar(-b,-mu,S,nu)
  out$mean = -out$mean
  return(out)
}


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

ftmvtmuvar = function(a,b,mu,S,nu){
  
  S = as.matrix(S)
  n = length(mu)
  s = sqrt(diag(S))
  a1 = a-mu
  b1 = b-mu
  as = a1/s
  bs = b1/s
  
  log2p = pmvt.genz(lower = a1,upper = b1,nu = nu,sigma = S,uselog2 = TRUE,N = 799)$Estimation
  p = 2^log2p
  
  if(n == 1){
    
    muY  = 0
    varY = 0
    kappa = log2ratio(s*gamma((nu+1)/2)/gamma(nu/2)/sqrt(nu*pi),log2p)
    
    av = (1+as^2/nu)^(-(nu-1)/2)
    bv = (1+bs^2/nu)^(-(nu-1)/2)
    
    if(nu == 1){
      muY = mu+s*log((1+bs^2)/(1+as^2))/(atan(bs)-atan(as))/2
    }else{
      muY = mu+kappa*nu/(nu-1)*(av-bv)
    }
    if (nu == 2){
      c = asinh(bs/sqrt(2))-asinh(as/sqrt(2))+(as*av-bs*bv)/sqrt(2)
      varY = log2ratio(c*S,log2p)-(muY-mu)^2
    }else{
      varY = ((nu-3)*mu*muY+(nu*S+mu^2)+nu*kappa*(a*av-b*bv))/(nu-2) - muY^2
    }
      return(list(mean = muY,EYY = varY + muY^2,varcov = varY))
  }
  
  qout = qfunT(a1,b1,S,nu)
  qa = qout$qa
  qb = qout$qb
  q = qa-qb
  muY = mu+S%*%log2ratio(q,log2p)
  
  if(n == 2){
    rho = S[1,2]/(s[1]*s[2])
    
    if(nu != 2){
      #  A faster way of computing nu/(nu-2)*mvtcdf(as/c,bs/c,R,nu-2), where c=sqrt(nu/(nu-2))
      p0 = (nu*p + sum(a1*qa) - sum(b1*qb))/(nu-2)
      
    }else{
      faux = function(w1){
        return(1/sqrt(2+w1^2)*
                 (atan((bs[2]-rho*w1)/sqrt(1-rho^2)/sqrt(2+w1^2)) -
                    atan((as[2]-rho*w1)/sqrt(1-rho^2)/sqrt(2+w1^2))))
      }
      p0 = integrate(faux,as[1],bs[1],abs.tol = 1e-10,rel.tol = 1e-10)$value/pi
    }
    Si = solve(S)
    c = c(t(a1)%*%Si%*%a1,
          t(c(a1[1],b1[2]))%*%Si%*%c(a1[1],b1[2]),
          t(c(b1[1],a1[2]))%*%Si%*%c(b1[1],a1[2]),
          t(b1)%*%Si%*%b1)
    if(nu==2){
      p1 = -log(2+c)
    }else{
      p1 = nu/(nu-2)*(1+c/nu)^(-nu/2+1)
    }
    
    h = sqrt(1-rho^2)*(p1[1]-p1[2]-p1[3]+p1[4])/(2*pi)
    W = matrix(0,2,2)
    
    #  W = mu*q'*S+S*(D-p0*eye(2))
    W[1,1] = S[1,1]*((a[1]+mu[1])*qa[1]-(b[1]+mu[1])*qb[1])+
      S[1,2]*((2*mu[1]+S[1,2]/S[2,2]*a1[2])*qa[2]-(2*mu[1]+S[1,2]/S[2,2]*b1[2])*qb[2])+S[1,1]*rho*h
    W[2,2] = S[2,2]*((a[2]+mu[2])*qa[2]-(b[2]+mu[2])*qb[2])+
      S[1,2]*((2*mu[2]+S[1,2]/S[1,1]*a1[1])*qa[1]-(2*mu[2]+S[1,2]/S[1,1]*b1[1])*qb[1])+S[2,2]*rho*h
    W[1,2] = S[1,1]*mu[2]*q[1]+S[2,2]*mu[1]*q[2]+
      S[1,2]*(a[1]*qa[1]+a[2]*qa[2]-b[1]*qb[1]-b[2]*qb[2])+s[1]*s[2]*h
    W[2,1] = W[1,2]
    
    varY = log2ratio(p0*S + W,log2p) + mu%*%t(mu) - muY%*%t(muY)
    
  }else{
    h = sqrt(nu)/(2*sqrt(pi)*gamma(nu/2))
    
    D =  lfun(a1,b1,nu*S,nu-2)*nu/(2*gamma(nu/2))*diag(n)
    
    for(i in 1:n){
      if(is.finite(a[i])){
        D[i,i] = D[i,i]+a[i]*qa[i]
      }
      if(is.finite(b[i])){
        D[i,i] = D[i,i]-b[i]*qb[i]
      }
    }
    
    for(i in 1:n){
      RR = S[-i,-i]-S[-i,i]%*%t(S[i,-i])/S[i,i]
      
      if(is.finite(a[i])){
        ma = mu[-i]+S[-i,i]/S[i,i]*a1[i]
        wa = hfun(a[-i],b[-i],ma,RR*(nu+as[i]^2),nu-1)*(1+as[i]^2/nu)^(-(nu-1)/2)*h/s[i]
      }else{
        wa = matrix(0,n-1,1)
      }
      if(is.finite(b[i])){
        mb = mu[-i]+S[-i,i]/S[i,i]*b1[i]
        wa = wa-hfun(a[-i],b[-i],mb,RR*(nu+bs[i]^2),nu-1)*(1+bs[i]^2/nu)^(-(nu-1)/2)*h/s[i]
      }
      
      D[i,-i] = wa
    }
    
    varY = log2ratio(S%*%(D-q%*%t(muY)),log2p)
  }
  
  varY = (varY + t(varY))/2
  
    return(list(mean = muY,EYY = varY + muY%*%t(muY),varcov = varY))
}

