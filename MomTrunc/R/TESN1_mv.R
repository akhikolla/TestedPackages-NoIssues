meanvarESNuni = function(a=-Inf,b=Inf,mu,Sigma,lambda,tau){
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
      Delta = s*lambda/sqrt(1+lambda^2)
      Omega  = cbind(rbind(Sigma,-Delta),rbind(-Delta,1))
      rownames(Omega) <- colnames(Omega)
      out = Kan.LRIC(a = c(a,-Inf),b = c(b,tautil),mu = c(mu,0),Sigma = Omega)
      return(list(mean = out$mean[-2],EYY = out$EYY[1,1],varcov = out$varcov[1,1]))
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
    #delta = lambda*s*eta
    F1N = (mu-mub)*F0N + Gamma*(dnorm(a,mu-mub,sqrt(Gamma)) - dnorm(b,mu-mub,sqrt(Gamma)))
    F2 = mu*F1 + Sigma*F0 + Sigma*(ifelse(a == -Inf,0,a*da) - ifelse(b == Inf,0,b*db)) + lambda*s*eta*F1N
    varY = c(F2/F0) - muY^2
    return(list(mean = muY,EYY = c(F2/F0), varcov = varY))
}
