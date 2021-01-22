####################################################################
#AUXILIAR CODES T
####################################################################

onlymeanT = function(a,b,mu,Sigma,nu)
{
  p = length(mu)
  nnu = nu/(nu-2)
  if(p==1){
    F0 = pent(b,mu,Sigma,nu) - pent(a,mu,Sigma,nu)
    nnusigma2 = nnu*Sigma
    ta = dent(a,mu,nnusigma2,nu-2)
    tb = dent(b,mu,nnusigma2,nu-2)
    F1 = mu*F0 + nnusigma2*(ta-tb)
    return(list(mean = as.numeric(F1/F0)))
  }
  
  bool1 = is.infinite(a)
  bool2 = is.infinite(b)
  
  if(sum(bool1*bool2) > 0){ #Does exist infinite pairs?
    
    if(sum(bool1*bool2) == p){ #All infinites?
      
      #NO TRUNCATION
      return(list(mean = mu))
      
    }else{
      
      return(withinfsT(a,b,mu,Sigma,nu))
      
    }
  }
  
  logF0 = pmvt.genz(lower = a-mu,upper = b-mu,nu = nu,sigma = Sigma,uselog2 = TRUE)$Estimation
  
  #Vectors ca and cb
  SSigma  = nnu*Sigma
  ssigma2 = diag(SSigma)
  ca = cb = rep(0,p)
  deltaA  = (nu - 2 + ((a - mu)^2)/diag(SSigma))/(nu - 1)
  deltaB  = (nu - 2 + ((b - mu)^2)/diag(SSigma))/(nu - 1)
  yA      = (a - mu)/diag(SSigma)
  yB      = (b - mu)/diag(SSigma)

  for(j in 1:p)
  {
    if(a[j]!=-Inf){
      ca[j] = log2prod(dent(a[j],mu[j],ssigma2[j],nu-2),
                       pmvt.genz(lower = a[-j] - (mu[-j] + yA[j]*SSigma[,j][-j]),
                                 upper = b[-j] - (mu[-j] + yA[j]*SSigma[,j][-j]),
                                 nu = nu - 1,
                                 sigma = deltaA[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),
                                 uselog2 = TRUE)$Estimation)
    }
    if(b[j]!= Inf){
      cb[j] = log2prod(dent(b[j],mu[j],ssigma2[j],nu-2),
                       pmvt.genz(lower = a[-j] - (mu[-j] + yB[j]*SSigma[,j][-j]),
                                 upper = b[-j] - (mu[-j] + yB[j]*SSigma[,j][-j]),
                                 nu = nu - 1,
                                 sigma = deltaB[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),
                                 uselog2 = TRUE)$Estimation)
    }
  }
  muY  = mu + log2ratio(SSigma%*%(ca - cb),logF0)
  return(list(mean = muY))
}


######################################################################################################################################
######################################################################################################################################

onlymeanT0 = function(a,b,mu,Sigma,nu)
{
  p = length(mu)
  nnu = nu/(nu-2)
  if(p==1){
    F0 = pent(b,mu,Sigma,nu) - pent(a,mu,Sigma,nu)
    nnusigma2 = nnu*Sigma
    ta = dent(a,mu,nnusigma2,nu-2)
    tb = dent(b,mu,nnusigma2,nu-2)
    F1 = mu*F0 + nnusigma2*(ta-tb)
    return(list(mean = as.numeric(F1/F0),logF00 = log2(F0)))
  }
  
  
  bool1 = is.infinite(a)
  bool2 = is.infinite(b)
  
  if(sum(bool1*bool2) > 0){ #Does exist infinite pairs?
    
    if(sum(bool1*bool2) == p){ #All infinites?
      
      #NO TRUNCATION
      return(list(mean = mu,logF00 = 0))
      
    }else{
      
      return(withinfs_meanT0(a,b,mu,Sigma,nu))
      
    }
  }

  logF0 = pmvt.genz(lower = a-mu,upper = b-mu,nu = nu,sigma = Sigma,uselog2 = TRUE)$Estimation
  
  #Vectors ca and cb
  SSigma  = nnu*Sigma
  ssigma2 = diag(SSigma)
  ca = cb = rep(0,p)
  deltaA  = (nu - 2 + ((a - mu)^2)/diag(SSigma))/(nu - 1)
  deltaB  = (nu - 2 + ((b - mu)^2)/diag(SSigma))/(nu - 1)
  yA      = (a - mu)/diag(SSigma)
  yB      = (b - mu)/diag(SSigma)

  for(j in 1:p)
  {
    if(a[j]!=-Inf){
      ca[j] = log2prod(dent(a[j],mu[j],ssigma2[j],nu-2),
                       pmvt.genz(lower = a[-j] - (mu[-j] + yA[j]*SSigma[,j][-j]),
                                 upper = b[-j] - (mu[-j] + yA[j]*SSigma[,j][-j]),
                                 nu = nu - 1,
                                 sigma = deltaA[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),
                                 uselog2 = TRUE)$Estimation)
    }
    if(b[j]!= Inf){
      cb[j] = log2prod(dent(b[j],mu[j],ssigma2[j],nu-2),
                       pmvt.genz(lower = a[-j] - (mu[-j] + yB[j]*SSigma[,j][-j]),
                                 upper = b[-j] - (mu[-j] + yB[j]*SSigma[,j][-j]),
                                 nu = nu - 1,
                                 sigma = deltaB[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),
                                 uselog2 = TRUE)$Estimation)
    }
  }
  muY  = mu + log2ratio(SSigma%*%(ca - cb),logF0)
  return(list(mean = muY,logF00 = logF0))
}
