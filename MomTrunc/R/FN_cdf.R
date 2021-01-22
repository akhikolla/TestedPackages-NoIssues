cdfFN = function(x,mu,Sigma)
{
  p = length(mu)
  if(p==1)
  {
    return(pnorm2(-x,x,mu,sqrt(Sigma)))
  }else
  {
    return(pmvnormt(lower = -x-mu,upper = x-mu,sigma = Sigma))
  }
}

