meanvarT16 = function(a,b,mu,Sigma,nu,omega = FALSE)
{
  
  p = length(mu)
  nnu = nu/(nu-2)
  if(p==1){
    if(nu==3){nu = 3.000000005}
    F0 = pent2(a,b,mu,Sigma,nu)
    nnusigma2 = nnu*Sigma
    ta = dent(a,mu,nnusigma2,nu-2)
    tb = dent(b,mu,nnusigma2,nu-2)
    L2 = pent2(a,b,mu,nnusigma2,nu-2)
    F1 = mu*F0 + nnusigma2*(ta-tb)
    F2 = mu*F1 + nnusigma2*(L2 + ifelse(a==-Inf,0,a*ta) - ifelse(b==Inf,0,b*tb))
    if(omega){
      return(list(mean = F1/F0,EYY = F2/F0,varcov = F2/F0 - (F1/F0)^2,omega = L2/F0*nu/(nu-2)))  
    }else{
      return(list(mean = F1/F0,EYY = F2/F0,varcov = F2/F0 - (F1/F0)^2))
    }
  }
  #GB = GenzBretz(maxpts = (p-1)*1e4, abseps = 1e-6, releps = 0)
  #print(GB$maxpts)
  #F0 = pmvt(lower = a-mu,upper = b-mu,df = nu,sigma = Sigma)[1]
  
  if(nu==3){nu = 3.01}
  
  logF0 = pmvt.genz(lower = a-mu,upper = b-mu,nu = nu,sigma = Sigma,uselog2 = TRUE,N = 799)$Estimation
  logF0nnu = pmvt.genz(lower = a-mu,upper = b-mu,nu = nu - 2,sigma = nnu*Sigma,N = 799,uselog2 = TRUE)$Estimation
  
  #Vectors ca and cb
  SSigma  = nnu*Sigma
  ssigma2 = diag(SSigma)
  ca = cb = rep(0,p)
  deltaA  = (nu - 2 + ((a - mu)^2)/diag(SSigma))/(nu - 1)
  deltaB  = (nu - 2 + ((b - mu)^2)/diag(SSigma))/(nu - 1)
  yA      = (a - mu)/diag(SSigma)
  yB      = (b - mu)/diag(SSigma)
  
  Daa = Dab = Dba = Dbb = matrix(0,p,p)
  #Paa = Pab = Pba = Pbb = matrix(0,p,p)
  
  Wa = Wb = matrix(0,p,p)
  
  nu0    = nu-1
  nnu0    = nu0/(nu0-2)
  
  seqq = seq_len(p)
  
  nf = genPxy(a,b,mu,Sigma,nu)
  
  for(j in seqq)
  {
    seqq_j = seq_len(p)[-j]
    #W matrix construction
    
    a0     = a[-j]
    b0     = b[-j]
    mu0a    = mu[-j] + yA[j]*SSigma[j,-j]
    Sigma0a = deltaA[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,-j])/ssigma2[j])
    
    mu0b    = mu[-j] + yB[j]*SSigma[j,-j]
    Sigma0b = deltaB[j]*(SSigma[-j,-j] - SSigma[j,-j]%*%t(SSigma[j,-j])/ssigma2[j])
    
    if(a[j]!=-Inf){
      
      for(i in 1:(p-1))
      {
        Daa[seqq_j[i],j] = dent(a0[i],mu0a[i],nnu0*diag(as.matrix(Sigma0a))[i],nu-3)
        Dab[seqq_j[i],j] = dent(b0[i],mu0a[i],nnu0*diag(as.matrix(Sigma0a))[i],nu-3)
      }
      
      logF00a  = pmvnormt(a[-j],b[-j],mu0a,Sigma0a,nu-1,uselog2 = TRUE)
      ca[j]    = log2prod(dent(a[j],mu[j],ssigma2[j],nu-2),logF00a)
      
      Wa[-j,j] = mu0a + log2ratio(nnu0*Sigma0a%*%(Daa[-j,j]*nf$Paa[-j,j] - Dab[-j,j]*nf$Pab[-j,j]),logF00a)
      Wa[j,j]  = a[j]
    }
    
    if(b[j]!= Inf){
      
      for(i in 1:(p-1))
      {
        Dba[seqq_j[i],j] = dent(a0[i],mu0b[i],nnu0*diag(as.matrix(Sigma0b))[i],nu-3)
        Dbb[seqq_j[i],j] = dent(b0[i],mu0b[i],nnu0*diag(as.matrix(Sigma0b))[i],nu-3)
      }
      
      logF00b  = pmvnormt(a[-j],b[-j],mu0b,Sigma0b,nu-1,uselog2 = TRUE)
      cb[j]    = log2prod(dent(b[j],mu[j],ssigma2[j],nu-2),logF00b)
      
      Wb[-j,j] = mu0b + log2ratio(nnu0*Sigma0b%*%(Dba[-j,j]*nf$Pba[-j,j] - Dbb[-j,j]*nf$Pbb[-j,j]),logF00b)
      Wb[j,j]  = b[j]
    }
    
  }
  
  muY  = mu + log2ratio(SSigma%*%(ca - cb),logF0)
  
  ratio0 = 2^(logF0nnu - logF0)
  Exx  = muY%*%t(mu) +  ratio0*SSigma +  log2ratio((Wa%*%diag(ca) - Wb%*%diag(cb))%*%SSigma,logF0)
  varY = Exx - muY%*%t(muY)
  if(omega){
    return(list(mean = muY,EYY = Exx,varcov = varY,omega = nu/(nu-2)*ratio0))
  }else{
    return(list(mean = muY,EYY = Exx,varcov = varY))
  }
}

# Upper -------------------------------------------------------------------

meanvarT16_upper = function(b,mu,Sigma,nu,omega = FALSE)
{
  if(nu==3){nu = 3.01}
  
  p = length(mu)
  nnu = nu/(nu-2)
  
  if(p==1){
    meanvarT16(-Inf,b,mu,Sigma,nu,omega)
  }
  #GB = GenzBretz(maxpts = (p-1)*1e4, abseps = 1e-6, releps = 0)
  #print(GB$maxpts)
  #F0 = pmvt(lower = a-mu,upper = b-mu,df = nu,sigma = Sigma)[1]
  
  logF0 = pmvt.genz(upper = b-mu,nu = nu,sigma = Sigma,uselog2 = TRUE,N = 799)$Estimation
  logF0nnu = pmvt.genz(upper = b-mu,nu = nu - 2,sigma = nnu*Sigma,uselog2 = TRUE,N = 799)$Estimation
  
  #Vectors ca and cb
  SSigma  = nnu*Sigma
  ssigma2 = diag(SSigma)
  cb = rep(0,p)
  deltaB  = (nu - 2 + ((b - mu)^2)/diag(SSigma))/(nu - 1)
  #yA      = (a - mu)/diag(SSigma)
  yB      = (b - mu)/ssigma2
  
  Dbb = matrix(0,p,p)
  #Paa = Pab = Pba = Pbb = matrix(0,p,p)
  
  Wb = matrix(0,p,p)
  
  nu0    = nu-1
  nnu0    = nu0/(nu0-2)
  
  seqq = seq_len(p)
  
  nf = genPxy_upper(b,mu,Sigma,nu)
  
  for(j in seqq)
  {
    seqq_j = seq_len(p)[-j]
    #W matrix construction
    
    b0     = b[-j]
    mu0b    = mu[-j] + yB[j]*SSigma[j,-j]
    Sigma0b = deltaB[j]*(SSigma[-j,-j] - SSigma[j,-j]%*%t(SSigma[j,-j])/ssigma2[j])
    
    if(b[j]!= Inf){
      
      for(i in 1:(p-1))
      {
        Dbb[seqq_j[i],j] = dent(b0[i],mu0b[i],nnu0*diag(as.matrix(Sigma0b))[i],nu-3)
      }
      
      logF00b  = pmvnormt(upper = b[-j],mean = mu0b,sigma = Sigma0b,nu = nu-1,uselog2 = TRUE)
      cb[j]    = log2prod(dent(b[j],mu[j],ssigma2[j],nu-2),logF00b)
      
      Wb[-j,j] = mu0b - log2ratio(nnu0*Sigma0b%*%(Dbb[-j,j]*nf$Pbb[-j,j]),logF00b)
      
      #Wb[-j,j] = mu0b + log2ratio(-nnu0*Sigma0b%*%(Dbb[-j,j]*nf$Pbb[-j,j]),logF00b)
      
      Wb[j,j]  = b[j]
    }
    
  }
  
  muY  = mu - log2ratio(SSigma%*%cb,logF0)
  ratio0 = 2^(logF0nnu - logF0)
  
  Exx  = muY%*%t(mu) +  ratio0*SSigma - log2ratio(Wb%*%diag(cb)%*%SSigma,logF0)
  
  #Exx = (Exx + t(Exx))/2
  varY = Exx - muY%*%t(muY)
  if(omega){
    return(list(mean = muY,EYY = Exx,varcov = varY,omega = nu/(nu-2)*ratio0))
  }else{
    return(list(mean = muY,EYY = Exx,varcov = varY))
  }
}

# Lower -------------------------------------------------------------------

meanvarT16_lower = function(a,mu,Sigma,nu,omega = FALSE)
{
  
  out = meanvarT16_upper(-a,-mu,Sigma,nu,omega)
  out$mean = -out$mean
  return(out)
  
}


# Finite ------------------------------------------------------------------

meanvarT16_finite = function(a,b,mu,Sigma,nu,omega = FALSE)
{
  if(nu==3){nu = 3.01}
  p = length(mu)
  nnu = nu/(nu-2)
  if(p==1){
    F0 = pent2(a,b,mu,Sigma,nu)
    nnusigma2 = nnu*Sigma
    ta = dent(a,mu,nnusigma2,nu-2)
    tb = dent(b,mu,nnusigma2,nu-2)
    F1 = mu*F0 + nnusigma2*(ta-tb)
    L2 = pent2(a,b,mu,nnusigma2,nu-2)
    F2 = mu*F1 + nnusigma2*(L2 + a*ta - b*tb)
    if(omega){
      return(list(mean = F1/F0,EYY = F2/F0,varcov = F2/F0 - (F1/F0)^2,omega = L2/F0*nu/(nu-2)))  
    }else{
    return(list(mean = F1/F0,EYY = F2/F0,varcov = F2/F0 - (F1/F0)^2))
      }
  }
  #GB = GenzBretz(maxpts = (p-1)*1e4, abseps = 1e-6, releps = 0)
  #print(GB$maxpts)
  #F0 = pmvt(lower = a-mu,upper = b-mu,df = nu,sigma = Sigma)[1]
  
  logF0 = pmvt.genz(lower = a-mu,upper = b-mu,nu = nu,sigma = Sigma,uselog2 = TRUE,N = 799)$Estimation
  logF0nnu = pmvt.genz(lower = a-mu,upper = b-mu,nu = nu - 2,sigma = nnu*Sigma,uselog2 = TRUE,N = 799)$Estimation
  
  #Vectors ca and cb
  SSigma  = nnu*Sigma
  ssigma2 = diag(SSigma)
  ca = cb = rep(0,p)
  deltaA  = (nu - 2 + ((a - mu)^2)/diag(SSigma))/(nu - 1)
  deltaB  = (nu - 2 + ((b - mu)^2)/diag(SSigma))/(nu - 1)
  yA      = (a - mu)/diag(SSigma)
  yB      = (b - mu)/diag(SSigma)
  
  Daa = Dab = Dba = Dbb = matrix(0,p,p)
  #Paa = Pab = Pba = Pbb = matrix(0,p,p)
  
  Wa = Wb = matrix(0,p,p)
  
  nu0    = nu-1
  nnu0    = nu0/(nu0-2)
  
  seqq = seq_len(p)
  
  nf = genPxy_finite(a,b,mu,Sigma,nu)
  
  for(j in seqq)
  {
    seqq_j = seq_len(p)[-j]
    #W matrix construction
    
    a0     = a[-j]
    b0     = b[-j]
    mu0a    = mu[-j] + yA[j]*SSigma[j,-j]
    Sigma0a = deltaA[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,-j])/ssigma2[j])
    
    mu0b    = mu[-j] + yB[j]*SSigma[j,-j]
    Sigma0b = deltaB[j]*(SSigma[-j,-j] - SSigma[j,-j]%*%t(SSigma[j,-j])/ssigma2[j])
      
      for(i in 1:(p-1))
      {
        Daa[seqq_j[i],j] = dent(a0[i],mu0a[i],nnu0*diag(as.matrix(Sigma0a))[i],nu-3)
        Dab[seqq_j[i],j] = dent(b0[i],mu0a[i],nnu0*diag(as.matrix(Sigma0a))[i],nu-3)
        Dba[seqq_j[i],j] = dent(a0[i],mu0b[i],nnu0*diag(as.matrix(Sigma0b))[i],nu-3)
        Dbb[seqq_j[i],j] = dent(b0[i],mu0b[i],nnu0*diag(as.matrix(Sigma0b))[i],nu-3)
      }
      
      logF00a  = pmvnormt(a[-j],b[-j],mu0a,Sigma0a,nu-1,uselog2 = TRUE)
      ca[j]    = log2prod(dent(a[j],mu[j],ssigma2[j],nu-2),logF00a)
      
      Wa[-j,j] = mu0a + log2ratio(nnu0*Sigma0a%*%(Daa[-j,j]*nf$Paa[-j,j] - Dab[-j,j]*nf$Pab[-j,j]),logF00a)
      Wa[j,j]  = a[j]
    

      logF00b  = pmvnormt(a[-j],b[-j],mu0b,Sigma0b,nu-1,uselog2 = TRUE)
      cb[j]    = log2prod(dent(b[j],mu[j],ssigma2[j],nu-2),logF00b)
      
      Wb[-j,j] = mu0b + log2ratio(nnu0*Sigma0b%*%(Dba[-j,j]*nf$Pba[-j,j] - Dbb[-j,j]*nf$Pbb[-j,j]),logF00b)
      Wb[j,j]  = b[j]
    
    
  }
  
  muY  = mu + log2ratio(SSigma%*%(ca - cb),logF0)
  #Exx  = muY%*%t(mu) +  log2ratio((2^logF0nnu*diag(p) + Wa%*%diag(ca) - Wb%*%diag(cb))%*%SSigma,logF0)
  #Exx  = muY%*%t(mu) +  log2ratio(2^logF0nnu*SSigma,logF0) +  log2ratio((Wa%*%diag(ca) - Wb%*%diag(cb))%*%SSigma,logF0)
  
  ratio0 = 2^(logF0nnu - logF0)
  
  Exx  = muY%*%t(mu) +  ratio0*SSigma +  log2ratio((Wa%*%diag(ca) - Wb%*%diag(cb))%*%SSigma,logF0)
  
  #Exx = (Exx + t(Exx))/2
  varY = Exx - muY%*%t(muY)
  if(omega){
    return(list(mean = muY,EYY = Exx,varcov = varY,omega = nu/(nu-2)*ratio0))
  }else{
  return(list(mean = muY,EYY = Exx,varcov = varY))
  }
}

# Test --------------------------------------------------------------------

# p = 5
# a  = seq(-0.9,0,length.out = p)
# b  = seq(0,0.9,length.out = p) + seq(0.1,0.6,length.out = p)
# mu = seq(0.25,0.75,length.out = p)*b
# s  = matrix(1.2*rnorm(p^2),p,p)
# Sigma = S = s%*%t(s)
# nu = 5
# 
# meanvarT16(a,b,mu,Sigma,5)
# 
# meanvarT16_finite(a,b,mu,Sigma,5)
# 
# a[c(1,3)] = -Inf
# meanvarT16(a,b,mu,Sigma,5)
# 
# b[c(3,5)] = Inf
# meanvarT16(a,b,mu,Sigma,5)
# meanvarT.Lin.LRIC(a,b,mu,Sigma,5)
# 
# a = rep(-Inf,p)
# 
# meanvarT16_upper(b,mu,Sigma,5)
# 
# meanvarT.Lin.RC(b,mu,Sigma,5)
# 
# meanvarT_upper(b,mu,Sigma,5)
# 
# 
# b = rep(Inf,p)
# 
# meanvarT16_lower(a,mu,Sigma,5)
# 
# meanvarT.Lin.LC(a,mu,Sigma,5)
# 
# meanvarT_lower(a,mu,Sigma,5)

# b = rep(Inf,p)
# 
# rbenchmark::benchmark(meanvarT16_lower(a,mu,Sigma,nu = 4),
#                       meanvarT_lower(a,mu,Sigma,nu=4),replications = 10)
# 
# p = 5
# a  = seq(-0.9,0,length.out = p)
# b  = seq(0,0.9,length.out = p) + seq(0.1,0.6,length.out = p)
# 
# rbenchmark::benchmark(meanvarT16_finite(a,b,mu,Sigma,nu = 5),
#                       meanvarT.Lin.IC(a,b,mu,Sigma,nu = 5),
#                       meanvarT_finite(a,b,mu,Sigma,nu=5),replications = 10)
