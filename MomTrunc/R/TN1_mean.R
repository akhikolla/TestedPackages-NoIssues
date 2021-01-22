onlymeanNuni = function(a=-Inf,b=Inf,mu,Sigma){
  if(a == -Inf & b==Inf){
    return(list(mean = mu))
  }
  s = sqrt(Sigma)
  F0    = pnorm2(a,b,mu,s)
  if(F0 < 1e-250){
    val = select(a,b,mu,s)
    return(list(mean = val))
  }
  da     = dnorm(a,mu,s)
  db     = dnorm(b,mu,s)
  F1     = mu*F0 + Sigma*(da - db)
  muY    = c(mu + Sigma*(da - db)/F0)
  return(list(mean = muY))
}
