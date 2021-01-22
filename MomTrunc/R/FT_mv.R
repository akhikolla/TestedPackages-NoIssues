#######################################################################
#FOLDED T
#######################################################################

meanvarFT = function(mu,Sigma,nu)
{
  n    = length(x = mu)
  muY  = matrix(data = NA,nrow = n,ncol = 1)
  varY = matrix(data = NA,nrow = n,ncol = n)
  if(nu<=1){
    return(list(muY = muY,varY = varY))
  }
  s   = sqrt(diag(Sigma))
  h   = mu/s
  c0  = sqrt(nu)*gamma((nu-1)/2)/(sqrt(pi)*gamma(nu/2))*(1+h^2/nu)^(-(nu-1)/2)
  c1  = 2*pt(q = h,df = nu)- 1
  muY = s*(c0 + h*c1)

  if(nu<=2){
    return(list(muY = muY,varY = varY))
  }
  R = Sigma/(s%*%t(s))

  for (i in 1:n){
    varY[i,i] = mu[i]^2 + s[i]^2*nu/(nu-2)
    for (j in seq_len(i-1)){
      r   = R[i,j]
      eta = sqrt(1 - r^2)
      Sigma   = matrix(data = c(1,r,r,1),nrow = 2,ncol = 2)
      pmt = pmvt.genz(lower=-Inf, upper=c(h[i],h[j]),nu=nu,sigma=Sigma)[[1]]
      p   = 4*pmt - c1[i] - c1[j] - 1
      zij = sqrt(nu-1)*(h[i]-r*h[j])/eta/sqrt(nu+h[j]^2)
      zji = sqrt(nu-1)*(h[j]-r*h[i])/eta/sqrt(nu+h[i]^2)
      cc  = (1+(h[i]^2+h[j]^2-2*r*h[i]*h[j])/eta^2/nu)^(-(nu-2)/2)
      r1  = r/(nu-2)
      v10 = pt(q = zij,df = nu-1)
      v20 = pt(q = zji,df = nu-1)
      varY[i,j] = s[i]*s[j]*(p*(h[i]*h[j]+nu*r1)+(h[i]-r1*h[j])*c0[j]*(2*v10-1)+
                               (h[j]-r1*h[i])*c0[i]*(2*v20-1)+2*nu/(nu-2)*eta/pi*cc)
      varY[j,i] = varY[i,j]
    }
  }
  varY2 = varY - muY%*%t(muY)
  return(list(muY = muY,EYY = varY,varY = varY2))
}
