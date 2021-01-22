ESN.NOTRUNC.onlymean = function(mu,Sigma,lambda,tau){
  p =  length(lambda)
  if(p==1){
    SL   = lambda^2
    mom1 = invmills(tau,0,sqrt(1+lambda^2))*lambda
    return(list(mean = mu + sqrt(Sigma)*mom1))
  }
  SL = sum(lambda^2)
  mom1 = invmills(tau,0,sqrt(1+SL))*lambda
  SS = sqrtm(Sigma)
  mom1 = SS%*%mom1 + mu
  return(list(mean = mom1))
}

# ESN.NOTRUNC = function(mu,Sigma,lambda,tau){
#   p = length(mu)
#   if(p==1){
#     s = sqrt(Sigma)
#     tautil = tau/sqrt(1+sum(lambda^2))
#     phi    = lambda/sqrt(1+sum(lambda^2))
#     eta    = invmills(tau,0,sqrt(1+lambda^2))
#     Gamma  = Sigma/(1+lambda^2)
#     mub    = lambda*tau*Gamma/s
#     F1     = mu + lambda*s*eta
#     F2 = mu*F1 + Sigma + lambda*s*eta*(mu-mub)
#     varY = F2 - F1^2
#     return(list(mean = round(F1,4),EYY = round(F2,4), varcov = round(varY,4)))
#   }
#   SS   = sqrtm(Sigma)
#   iSS  = solve(SS)
#   #aux
#   tautil = tau/sqrt(1+sum(lambda^2))     #ok #dif
#   Phi    = diag(p) + lambda%*%t(lambda)  #ok #dif
#   Gamma  = SS%*%solve(Phi)%*%SS
#   Gamma = (t(Gamma) + Gamma)/2
#   iGamma = iSS%*%Phi%*%iSS
#   varphi = iSS%*%lambda            #ok
#   mub    = tau*Gamma%*%varphi
#   eta = invmills(tau,0,sqrt(1+sum(lambda^2)))
#   delta  = eta*Sigma%*%varphi                    #ok
#
#   #ESN part
#   muY = mu + delta
#   varY = matrix(NA,p,p)
#   EXX = mu%*%t(muY) + delta%*%t(mu - mub) + Sigma
#   EXX = (EXX + t(EXX))/2
#   varY = EXX - muY%*%t(muY)
#   return(list(mean = round(muY,4),EYY = round(EXX,4), varcov = round(varY,4)))
# }

#
#
# library(mvtnorm)
# meanvarESN.infs(mu,Sigma,lambda,tau)
#
# #MGF.noinfs(a = -10^6,b=10^6,mu,Sigma,lambda,tau)
#
# meanvarESN3(mu = mu,Sigma = Sigma,lambda = lambda,tau = tau)
#
# library(microbenchmark)
# library(ggplot2)
# compare <- microbenchmark(meanvarESN.infs(mu,Sigma,lambda,tau),
#                           meanvarESN3(mu = mu,Sigma = Sigma,lambda = lambda,tau = tau),
#                           times = 50)
# autoplot(compare)
#
#
# meanvarESN.infs(2,5,-1,2)
# meanvarESN3(mu = 2,Sigma = 5,lambda = -1,tau = 2)
#
# meanvarESNuni(a=0,mu = 2,Sigma = 5,lambda = -1,tau = 2)


#meanvarESNuni(a = NULL,b = NULL,mu = 2,Sigma = 5,lambda = -1,tau = 2)
