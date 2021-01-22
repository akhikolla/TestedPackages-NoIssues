onlymeanESNuni = function(a=-Inf,b=Inf,mu,Sigma,lambda,tau){
  s = sqrt(Sigma)
  tautil = tau/sqrt(1+sum(lambda^2))
  if(tautil< -35){
    #print("normal aproximation")
    Gamma  = Sigma/(1+lambda^2)
    mub    = lambda*tau*Gamma/s
    return(meanvarNuni(a,b,mu-mub,Gamma))
  }
  F0     = AcumESN(b,mu,Sigma,lambda,tau) - AcumESN(a,mu,Sigma,lambda,tau)
  if(F0 < 1e-250){
    val = select(a,b,mu,s)
    return(list(mean = val))
  }
  phi    = lambda/sqrt(1+sum(lambda^2))
  eta    = invmills(tau,0,sqrt(1+lambda^2))
  Gamma  = Sigma/(1+lambda^2)
  mub    = lambda*tau*Gamma/s
  F0N    = pnorm2(a,b,mu-mub,sqrt(Gamma))
  da     = dmvESN1(a,mu,Sigma,lambda,tau)
  db     = dmvESN1(b,mu,Sigma,lambda,tau)
  F1     = mu*F0 + Sigma*(da - db) + lambda*s*eta*F0N
  muY = c(F1/F0)
  return(list(mean = muY))
}
