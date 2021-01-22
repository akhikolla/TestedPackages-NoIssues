
#######################################################################
#FOLDED NORMAL
#######################################################################

meanvarFN = function(mu,Sigma)
{
  n    = length(x = mu)
  muY  = matrix(data = NA,nrow = n,ncol = 1)
  varY = matrix(data = NA,nrow = n,ncol = n)
  s   = sqrt(diag(Sigma))
  h   = mu/s
  c0  = dnorm(h)
  c1  = 2*pnorm(h)- 1
  muY = s*(2*c0 + h*c1)
  R = Sigma/(s%*%t(s))

  h1 = h%*%matrix(1,1,n)
  A  = (h1 - R*t(h1))/sqrt(2*(1-R^2))
  diag(A) = 0
  gam = (2*pnorm(sqrt(2)*A) - 1)*h%*%t(c0)
  for (i in 1:n){
    varY[i,i] = mu[i]^2 + s[i]^2
    for (j in seq_len(i-1)){
      r   = R[i,j]
      eta = sqrt(1 - r^2)
      Sigma   = matrix(data = c(1,r,r,1),nrow = 2,ncol = 2)
      pmt = pmvnorm(lower=-Inf, upper=c(h[i],h[j]),corr=Sigma,algorithm = TVPACK)[1]
      p   = 4*pmt - c1[i] - c1[j] - 1
      c = sqrt(h[i]^2 + h[j]^2 - 2*r*h[i]*h[j])/eta
      varY[i,j] = s[i]*s[j]*(p*(h[i]*h[j]+r)+2*gam[i,j]+2*gam[j,i]+4*eta/sqrt(2*pi)*dnorm(c))
      varY[j,i] = varY[i,j]
    }
  }
  varY2 = varY - muY%*%t(muY)
  return(list(muY = muY,EYY = varY,varY = varY2))
}
