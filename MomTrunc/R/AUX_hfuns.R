hfunrec2 = function(k,a,b,mu,S,nu){
  n = length(mu)
  if(n==1){
    y = hfun1(k,a,b,mu,S,nu)
    return(list(index = 0:k,y=y))
  }
  # %
  # %   Create a matrix of binomial coefficients of C(i,j) = nchoosek(i-1,j-1)
  # %
  C = matrix(0,n+k+1,n+1)
  C[,1] = 1
  for(j in 2:(n+1)){
    C[j:nrow(C),j] = cumsum(C[(j-1):(n+k),j-1])
  }
  y = rep(NaN,C[nrow(C),ncol(C)])
  s = sqrt(diag(S))
  as = (a-mu)/s
  bs = (b-mu)/s
  y[1] = hfun0(a,b,mu,S,nu)
  
  index = matrix(0,1,n)
  
  if(k>0){
    infind = is.infinite(a) | is.infinite(b)  # pair with infinity elements
    nu1 = nu+n-sum(infind)
    # %
    # %   G is a matrix with n columms.  Column i of G is used to store G_{kappa}^i for 
    # %   different values of kappa.  Each column has nchoosek(n+k-2,n-1) number of elements.
    # %
    ca = (1+as^2)^(-(nu+n-2)/2)/(nu+n-2)
    cb = (1+bs^2)^(-(nu+n-2)/2)/(nu+n-2)
    glen = C[n+k-1,n]  #% nchoosek(n+k-2,n-1)
    Ga = Gb = matrix(0,glen,n)
    
    
    for(i in 1:n){
      v = S[-i,-i]-S[-i,i]%*%t(S[-i,i])/S[i,i]
      ma = mu[-i]+S[-i,i]*as[i]/s[i]
      mb = mu[-i]+S[-i,i]*bs[i]/s[i]
      sa = (1+as[i]^2)*v
      sb = (1+bs[i]^2)*v
      if(is.finite(a[i])){
        Ga[,i] = ca[i]*hfunrec20(k-1,a[-i],b[-i],ma,sa,nu-1)
        #Ga[,i] = ca[i]*rep(1,(n-1+k-1+1)*(n-1+1))
      }
      if(is.finite(b[i])){
        Gb[,i] = cb[i]*hfunrec20(k-1,a[-i],b[-i],mb,sb,nu-1)
        #Gb[,i] = cb[i]*rep(1,(n-1+k-1+1)*(n-1+1))
      }
    }
    y[seq(n+1,2,by = -1)] = mu*y[1]+S%*%(Ga[1,]-Gb[1,])   #% H_{kappa} with |kappa|=1
    
    
    E = diag(n)
    start0 = 0         # Start index of |kappa|-2  
    start1 = 1         # Start index of |kappa|-1
    ind0 = matrix(0,1,n)
    ind1 = E[,rev(1:n)]
    index = rbind(index,ind1)
    
    if(k>1){
      Q  = rep(0,C[n+k-2,n])
      c  = solve(S,mu)
      cc = c(1 + t(mu)%*%c)
      
      for(ii in 2:k){       # |kappa|=ii
        #  Update Q  
        for(j in 1:nrow(ind0)){
          i1 = start1 + colexind(matrix(1,n,1)%*%ind0[j,]+E,C)
          Q[j] = cc*y[start0+j]-t(c)%*%y[i1]
          for(i in 1:n){
            i1 = C[n+ii-3-ind0[j,i],n]+colexind(ind0[j,-i],C)
            if(is.finite(a[i])){
              Q[j] = Q[j]+a[i]^(ind0[j,i]+1)*Ga[i1,i]
            }
            if(is.finite(b[i])){
              Q[j] = Q[j]-b[i]^(ind0[j,i]+1)*Gb[i1,i]
            }
          }
          Q[j] = Q[j]/(nu-ii)
        }
        ind0 = ind1
        start0 = start1
        start1 = start1 + nrow(ind1)
        ind1 = colex(ii,n)
        index = rbind(index,ind1)
        #ind1 = test2(ii,n)
        
        #colex?????
        #%  Update H_{kappa}
        
        
        for(j in 1:nrow(ind1)){
          if(sum(ind1[j,infind])<nu1){
            i = min(which(ind1[j,]>0))  # first nonzero index
            kappa = ind1[j,] - E[i,]
            i1 = start0 + colexind(kappa,C)
            y[start1+j] = mu[i]*y[i1]
            for(l in 1:n){
              i2 = C[n+ii-kappa[l]-2,n]+colexind(kappa[-l],C)
              if(is.finite(a[l])){
                y[start1+j] = y[start1+j]+S[i,l]*a[l]^kappa[l]*Ga[i2,l]
              }
              if(is.finite(b[l])){
                y[start1+j] = y[start1+j]-S[i,l]*b[l]^kappa[l]*Gb[i2,l]
              }
              if(kappa[l]>0)
                y[start1+j] = y[start1+j]+S[i,l]*kappa[l]*Q[colexind(kappa-E[l,],C)]
            }
          }
        }
      }
    }
  }
  return(list(index = index,y = y))
}

######################################################################################################%

hfunrec20 = function(k,a,b,mu,S,nu){
  n = length(mu)
  if(n==1){
    y = hfun1(k,a,b,mu,S,nu)
    return(y)
  }
  # %
  # %   Create a matrix of binomial coefficients of C(i,j) = nchoosek(i-1,j-1)
  # %
  C = matrix(0,n+k+1,n+1)
  C[,1] = 1
  for(j in 2:(n+1)){
    C[j:nrow(C),j] = cumsum(C[(j-1):(n+k),j-1])
  }
  y = rep(NaN,C[nrow(C),ncol(C)])
  s = sqrt(diag(S))
  as = (a-mu)/s
  bs = (b-mu)/s
  y[1] = hfun0(a,b,mu,S,nu)
  
  if(k>0){
    infind = is.infinite(a) | is.infinite(b)  # pair with infinity elements
    nu1 = nu+n-sum(infind)
    # %
    # %   G is a matrix with n columms.  Column i of G is used to store G_{kappa}^i for 
    # %   different values of kappa.  Each column has nchoosek(n+k-2,n-1) number of elements.
    # %
    ca = (1+as^2)^(-(nu+n-2)/2)/(nu+n-2)
    cb = (1+bs^2)^(-(nu+n-2)/2)/(nu+n-2)
    glen = C[n+k-1,n]  #% nchoosek(n+k-2,n-1)
    Ga = Gb = matrix(0,glen,n)
    
    
    for(i in 1:n){
      v = S[-i,-i]-S[-i,i]%*%t(S[-i,i])/S[i,i]
      ma = mu[-i]+S[-i,i]*as[i]/s[i]
      mb = mu[-i]+S[-i,i]*bs[i]/s[i]
      sa = (1+as[i]^2)*v
      sb = (1+bs[i]^2)*v
      if(is.finite(a[i])){
        Ga[,i] = ca[i]*hfunrec20(k-1,a[-i],b[-i],ma,sa,nu-1)
        #Ga[,i] = ca[i]*rep(1,(n-1+k-1+1)*(n-1+1))
      }
      if(is.finite(b[i])){
        Gb[,i] = cb[i]*hfunrec20(k-1,a[-i],b[-i],mb,sb,nu-1)
        #Gb[,i] = cb[i]*rep(1,(n-1+k-1+1)*(n-1+1))
      }
    }
    y[seq(n+1,2,by = -1)] = mu*y[1]+S%*%(Ga[1,]-Gb[1,])   #% H_{kappa} with |kappa|=1
    
    if(k>1){
      E = diag(n)
      start0 = 0         # Start index of |kappa|-2  
      start1 = 1         # Start index of |kappa|-1
      ind0 = matrix(0,1,n)
      ind1 = E[,rev(1:n)]
      Q  = rep(0,C[n+k-2,n])
      c  = solve(S,mu)
      cc = c(1 + t(mu)%*%c)
      
      for(ii in 2:k){       # |kappa|=ii
        #  Update Q  
        for(j in 1:nrow(ind0)){
          i1 = start1 + colexind(matrix(1,n,1)%*%ind0[j,]+E,C)
          Q[j] = cc*y[start0+j]-t(c)%*%y[i1]
          for(i in 1:n){
            i1 = C[n+ii-3-ind0[j,i],n]+colexind(ind0[j,-i],C)
            if(is.finite(a[i])){
              Q[j] = Q[j]+a[i]^(ind0[j,i]+1)*Ga[i1,i]
            }
            if(is.finite(b[i])){
              Q[j] = Q[j]-b[i]^(ind0[j,i]+1)*Gb[i1,i]
            }
          }
          Q[j] = Q[j]/(nu-ii)
        }
        ind0 = ind1
        start0 = start1
        start1 = start1 + nrow(ind1)
        ind1 = colex(ii,n)
        #ind1 = test2(ii,n)
        
        #colex?????
        #%  Update H_{kappa}
        
      
        for(j in 1:nrow(ind1)){
          if(sum(ind1[j,infind])<nu1){
            i = min(which(ind1[j,]>0))  # first nonzero index
            kappa = ind1[j,] - E[i,]
            i1 = start0 + colexind(kappa,C)
            y[start1+j] = mu[i]*y[i1]
            for(l in 1:n){
              i2 = C[n+ii-kappa[l]-2,n]+colexind(kappa[-l],C)
              if(is.finite(a[l])){
                y[start1+j] = y[start1+j]+S[i,l]*a[l]^kappa[l]*Ga[i2,l]
              }
              if(is.finite(b[l])){
                y[start1+j] = y[start1+j]-S[i,l]*b[l]^kappa[l]*Gb[i2,l]
              }
              if(kappa[l]>0)
                y[start1+j] = y[start1+j]+S[i,l]*kappa[l]*Q[colexind(kappa-E[l,],C)]
            }
          }
        }
      }
    }
  }
  return(y)
}

######################################################################################################%

colexind = function(kk,C){
  if(!is.matrix(kk)){kk = matrix(kk,ncol = length(kk))}
  m = nrow(kk)
  nn = ncol(kk)
  y1 = rep(1,m)
  for(jj in 1:m){
    s = cumsum(kk[jj,seq(nn,1,by = -1)])
    for(r in seq_len(nn-1)){
      if(kk[jj,r]>0){
        n1 = nn-r+1
        y1[jj] = y1[jj] + C[s[n1]+n1,n1] - C[s[n1-1]+n1,n1]
      }
    }
  }
  return(y1)
}

######################################################################################################%

# %
# %   This program computes H_k = int_a^b x^k[1+(x-mu)^2/v]^{-(nu+1)/2}dx
# %

hfun1 = function(k,a,b,mu,v,nu){
  
  fab = is.finite(a) && is.finite(b)
  s = sqrt(v)
  as = (a-mu)/s
  bs = (b-mu)/s
  av = 1/(1+as^2)^((nu-1)/2)
  bv = 1/(1+bs^2)^((nu-1)/2)
  y = matrix(NaN,k+1,1)
  # %
  # %  Compute H_0
  # %  
  if(nu==0 && fab){
    y[1] = s*(asinh(bs)-asinh(as))
  }else{
    if(nu<0 && fab){
      y[1] = ((nu+1)*hfun1(0,a,b,mu,v,nu+2)+(a-mu)/(1+as^2)^((nu+1)/2)-(b-mu)/(1+bs^2)^((nu+1)/2))/nu
    }else{
      if(nu>0){
        y[1] = s*beta(1/2,nu/2)*(pt(bs*sqrt(nu),nu)-pt(as*sqrt(nu),nu))
      }
    }
  }
  if(k>=1){
    if(fab){
      if(nu==1){
        y[2] = mu*y[1]+v/2*log((1+bs^2)/(1+as^2))
      }else{
        y[2] = mu*y[1]+v/(nu-1)*(av-bv)
      }
      for(i in seq_len(k-1)){
        if(i!=nu-1){
          y[i+2] = ((nu-1-2*i)*mu*y[i+1]+i*(mu^2+v)*y[i]+v*(a^i*av-b^i*bv))/(nu-1-i)
        }else{
          y1 = hfun1(i-1,a,b,mu,v,i-1)
          y[i+2] = y[i+1]*mu+v*y1[length(y1)]+(a^i/(1+as^2)^(i/2)-b^i/(1+bs^2)^(i/2))*v/i
        }             
      }
    }else{
      if(nu>1){
        y[2] = mu*y[1]+v/(nu-1)*(av-bv)
        for(i in seq_len(k-1)){
          if(i<nu-1){
            y[i+2] = (nu-1-2*i)*mu*y[i+1]+i*(mu^2+v)*y[i]
            if(is.finite(a)){
              y[i+2] = y[i+2]+v*a^i*av
            }
            if(is.finite(b)){
              y[i+2] = y[i+2]-v*b^i*bv
            }
            y[i+2] = y[i+2]/(nu-1-i)
          }
        }
      }
    }
  }
  return(y)
}

######################################################################################################%

# %
# %   This program computes H_{0_n}(a,b;mu,S,nu) = \int_a^b [1+(x-mu)'S^{-1}(x-mu]^{-(nu+n)/2}dx.
# %

hfun0 = function(a,b,mu,S,nu){
  n = length(mu)
  n1 = sum(is.infinite(a) & is.infinite(b))  # number of infinite pairs
  s = sqrt(diag(S))
  R = S/(s%*%t(s))
  as = (a-mu)/s
  bs = (b-mu)/s
  
  if(nu==0){
    
    if(n1==0){
      y = NaN
    }else{
      if(n==1){
        y = s*(asinh(bs)-asinh(as))
      }else{
        if(n==2){
          
          if(is.infinite(a[1]) && is.infinite(b[1])){
            y = pi*sqrt(det(S))*(asinh(bs[2])-asinh(as[2]))
          }else{
            
            if(is.infinite(a[2]) && is.infinite(b[2])){
              y = pi*sqrt(det(S))*(asinh(bs[1])-asinh(as[1]))
            }else{
              r = R[1,2]
              faux = function(w){
                return((atan((bs[2]-r*w)/sqrt((1-r^2)*(1+w^2)))-atan((as[2]-r*w)/sqrt((1-r^2)*(1+w^2))))/sqrt(1+w^2))
              }
              y = sqrt(det(S))*integrate(faux,as[1],bs[1],abs.tol = 1e-12,rel.tol = 1e-12)$value
            }
          }
        }else{
          ii = which(is.finite(a) & is.finite(b),TRUE) #find
          y = htmuk(ii[1],0,a,b,mu,S,0)
        }
      }
    }
  }else{
    
    if(nu>0){
      c = sqrt(det(S))*pi^(n/2)*gamma(nu/2)/gamma((nu+n)/2)
      y = pmvt.genz(lower = sqrt(nu)*as,upper = sqrt(nu)*bs,nu = nu,sigma = R)[[1]]
      y = c*y
    }else{
      # nu<0
      y = (nu+n)/nu*hfun0(a,b,mu,S,nu+2)
      if(n==1){
        if(is.finite(a)){
          y = y+(a-mu)/(1+as^2)^((nu+n)/2)/nu
        }
        if(is.finite(b)){
          y = y-(b-mu)/(1+bs^2)^((nu+n)/2)/nu
        }
      }else{
        for(j in 1:n){
          mua = mu+S[,j]*as[j]/s[j]
          mub = mu+S[,j]*bs[j]/s[j]
          Sj = S[-j,-j] -S[,-j]*S[-j,]/S[j,j]
          if(is.finite(a[j])){
            y = y+(a[j]-mu[j])/(1+as[j]^2)^((nu+n)/2)/nu*hfun0(a[-j],b[-j],mua[-j],(1+as[j]^2)*Sj,nu+1)
          }
          if(is.finite(b[j])){
            y = y-(b[j]-mu[j])/(1+bs[j]^2)^((nu+n)/2)/nu*hfun0(a[-j],b[-j],mub[-j],(1+bs[j]^2)*Sj,nu+1)
          }
        }
      }
    }
  }
  return(y)
}

######################################################################################################%

# %
# %   htmuk.m   Date: 5/10/2018
# %   This Matlab function computes H_{ke_i}(a,b;mu,S,nu)
# %   using a direct integration approach.  This function is needed
# %   because the recurrence approach fails when k=nu.
# %
htmuk = function(i,k,a,b,mu,S,nu){
  
  n = length(mu)
  cc = pi^((n-1)/2)*gamma((nu+1)/2)*sqrt(det(S))/gamma((nu+n)/2)
  s = sqrt(diag(S))
  as = (a-mu)/s
  bs = (b-mu)/s
  ind = 1:n
  ind = ind[-i]
  infpair = is.infinite(a[-i]) & is.infinite(b[-i])
  ind = ind[-infpair]
  a1  = as[ind]
  b1  = bs[ind]
  S   = S/(s%*%t(s))
  c   = S[i,ind]
  S0 = S[ind,ind] - t(c)%*%c
  s0 = sqrt(diag(S0))
  R = S0/(s0%*%t(s0))
  n1 = length(ind)
  integd = function(x){
    y1 = (1+x^2)^(-(nu+1)/2)
    if(n1>0){
      mut = c*x
      d = sqrt((1+x^2)/(nu+1))
      y1 = y1*pmvt(lower = (a1-mut)/(d*s0),upper = (b1-mut)/(d*s0),df = nu+1,sigma = R)[1]
      #y1 = y1*0.5
    }
  }
  faux1 = function(xp){
    return(sapply(xp,integd)*(mu[i]+s[i]*xp)^k)
  }
  y = cc*integrate(faux1,as[i],bs[i],abs.tol = 1e-10,rel.tol = 1e-10)$value
  #y = cc*integrate(faux1,as[i],bs[i])$value
  return(y)
}

# 
# htmuk(i,k,a,b,mu,S,nu)
