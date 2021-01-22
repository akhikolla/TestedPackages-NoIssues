cdfFESN = function(x,mu,Sigma,lambda,tau)
{
  return(pmvESN(lower = -x-mu,upper = x-mu,Sigma = Sigma,lambda = lambda,tau = tau))
}
