cdfFT = function(x,mu,Sigma,nu)
{
  p = length(mu)
  if(p==1)
  {
    return(pent2(-x,x,mu,Sigma,nu))
  }else
  {
    return(pmvt.genz(lower = -x-mu,upper = x-mu,nu = nu,sigma = Sigma)$Estimation)
  }
}

