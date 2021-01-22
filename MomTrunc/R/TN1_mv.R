meanvarNuni = function(a=-Inf,b=Inf,mu,Sigma){
  if(a == -Inf & b==Inf){
    return(list(mean = mu,EYY = Sigma + mu^2, varcov = Sigma))
  }
  s = sqrt(Sigma)
  F0    = pnorm2(a,b,mu,s)
  if(F0 < 1e-250){
    val = select(a,b,mu,s)
    return(list(mean = val,EYY = 0.0001 + val^2, varcov = 0.0001))
    }
  da     = dnorm(a,mu,s)
  db     = dnorm(b,mu,s)
  F1     = mu*F0 + Sigma*(da - db)
  muY    = c(mu + Sigma*(da - db)/F0)
  if(a == -Inf){a = 0}
  if(b ==  Inf){b = 0}
  EYY = mu*muY + Sigma + Sigma*(a*da - b*db)/F0
  varY = EYY - muY^2
  return(list(mean = muY,EYY = EYY, varcov = varY))
}
