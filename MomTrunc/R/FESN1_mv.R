# FOLDED UNIVARIATE ESN
meanvarFESN_uni = function(mu,Sigma,lambda,tau){
  s = sqrt(Sigma)
  slam   = sqrt(1+lambda^2)
  tautil = tau/slam
  if(tautil< -35){
    #print("normal aproximation")
    Gamma  = Sigma/(1+lambda^2)
    mub    = lambda*tau*Gamma/s
    return(meanvarFN(mu-mub,Gamma))
  }
  phi    = lambda/slam
  eta = invmills(tau,0,slam)
  gamma = sqrt(Sigma/slam^2)
  mub = lambda*tau*s/slam^2
  m = mu - mub
  FFF = AcumESN(0,mu,Sigma,lambda,tau)
  eee = dnorm(0,m,gamma)
  FNNN = pnorm(0,m,gamma)
  ddd = dmvESN1(0,mu,Sigma,lambda,tau)
  first = mu*(1-2*FFF) + 2*Sigma*ddd + lambda*eta*s*(1-2*FNNN)
  second = mu^2 + Sigma + eta*lambda*s*(m + mu)
  third = (1-2*FFF)*mu*(mu^2 + 3*Sigma) + 2*ddd*s*(s*mu^2 + 2*s^3) + eta*lambda*s*(2*eee*gamma^2*(m + mu) + (1-2*FNNN)*(m^2 + gamma^2 + mu*(m + mu) + 2*Sigma))
  fourth = mu^4 + 6*mu^2*Sigma + 3*Sigma^2 + eta*lambda*s*(m^3 + m^2*mu + m*(3*gamma^2 + mu^2 + 3*Sigma) + mu*(gamma^2 + mu^2 + 5*Sigma))
  return(list(muY=first,EYY=second,varcov=second-first^2))
}

onlymeanFESN_uni = function(mu,Sigma,lambda,tau){
  s = sqrt(Sigma)
  slam   = sqrt(1+lambda^2)
  tautil = tau/slam
  phi    = lambda/slam
  eta = invmills(tau,0,slam)
  gamma = sqrt(Sigma/slam^2)
  mub = lambda*tau*s/slam^2
  m = mu - mub
  FFF = AcumESN(0,mu,Sigma,lambda,tau)
  eee = dnorm(0,m,gamma)
  FNNN = pnorm(0,m,gamma)
  ddd = dmvESN1(0,mu,Sigma,lambda,tau)
  first = mu*(1-2*FFF) + 2*Sigma*ddd + lambda*eta*s*(1-2*FNNN)
  return(first)
}


# mu=3; Sigma=5; lambda = 3; tau=2
# a1 = rESN(n = 50000000,mu,Sigma,lambda,tau)
# b1 = abs(a1)
# r1 = mean(b1)
# r2 = mean(b1^2)
# r3 = mean(b1^3)
# r4 = mean(b1^4)
#
# r1;r2;r3;r4
# first;second;third;fourth
#
# # #################
# meanvarFESN_uni(mu,Sigma,lambda,tau)
