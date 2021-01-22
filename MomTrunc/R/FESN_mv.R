meanvarFESNopt = function(mu,Sigma,lambda,tau)
{
  p = length(mu)
  if(p==1){
    return(meanvarFESN_uni(mu,Sigma,lambda,tau))
  }
  tautil = tau/sqrt(1+sum(lambda^2))
  if(tautil < -35){
    #print("normal aproximation")
    Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
    return(meanvarFN(mu - tautil*Delta,Sigma = Sigma - Delta%*%t(Delta)))
  }
  muY  = matrix(data = NA,nrow = p,ncol = 1)
  EYY = matrix(data = NA,nrow = p,ncol = p)
  varY = matrix(data = NA,nrow = p,ncol = p)

  varphi = sqrtm(solve(Sigma))%*%lambda
  for (i in 1:p){
    tilSj     = Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i])
    cj        = (1 + t(varphi[-i])%*%tilSj%*%varphi[-i])^(-0.5)
    tilvarphi = varphi[i] + (Sigma[-i,i]/Sigma[i,i])%*%varphi[-i]
    out1      = meanvarFESN_uni(mu = mu[i],Sigma = Sigma[i,i],lambda = sqrt(Sigma[i,i])*cj*tilvarphi,tau = cj*tau)
    muY[i]    = out1[[1]]
    EYY[i,i]  = out1[[2]]
    varY[i,i] = out1[[3]]
    for (j in seq_len(i-1)){
      #Here we compute the new parameters for the bivariate partition Y = (Yi,Yj)
      mu.ij      = mu[c(i,j)]
      Sigma.ij   = Sigma[c(i,j),c(i,j)]
      iSigma.ij  = solve(Sigma.ij)
      SS.ij  = sqrtm(Sigma.ij)
      Sigma22.1  = Sigma[-c(i,j),-c(i,j)] - Sigma[-c(i,j),c(i,j)]%*%iSigma.ij%*%Sigma[c(i,j),-c(i,j)]
      c12        = as.numeric((1 + t(varphi[-c(i,j)])%*%Sigma22.1%*%varphi[-c(i,j)])^(-0.5))
      tilvarphi1 = varphi[c(i,j)] + iSigma.ij%*%Sigma[c(i,j),-c(i,j)]%*%varphi[-c(i,j)]
      lambda.ij  = c12*SS.ij%*%tilvarphi1
      tau.ij     = c12*tau
      EYY[i,j] = EYY[j,i] = Exixj(mu.ij,Sigma.ij,lambda.ij,tau.ij)
    }
  }
  varY = EYY - muY%*%t(muY)
  return(list(muY = muY,EYY = EYY,varY = varY))
}


Exixj = function(mu,Sigma,lambda,tau)
{
  SS = sqrtm(Sigma)
  iSS = solve(SS)
  Phi        = diag(2) + lambda%*%t(lambda)
  Gamma      = SS%*%solve(Phi)%*%SS
  Gammaneg   = diag(c(-1,1))%*%Gamma%*%diag(c(-1,1))
  Sigmaneg   = diag(c(-1,1))%*%Sigma%*%diag(c(-1,1))
  tautil     = tau/sqrt(1+sum(lambda^2))
  iGamma     = iSS%*%Phi%*%iSS
  varphi  = iSS%*%lambda            #ok
  mub        = tau*Gamma%*%varphi
  eta = invmills(tau,0,sqrt(1+sum(lambda^2)))

  delta  = eta*Sigma%*%varphi
  m = mu - mub

  PhiESN1 = pmvESN(upper = c(0,0),mu = c(-mu[1],mu[2]),Sigma = Sigmaneg,lambda = c(-lambda[1],lambda[2]),tau = tau)[1]
  PhiESN2 = pmvESN(upper = c(0,0),mu = c(mu[1],-mu[2]),Sigma = Sigmaneg,lambda = c(lambda[1],-lambda[2]),tau = tau)[1]
  rownames(Gammaneg) <- colnames(Gammaneg)
  Phi1 = pmvn.genz(upper = c(0,0),mean = c(-m[1],m[2]),sigma = Gammaneg)$Estimation
  Phi2 = pmvn.genz(upper = c(0,0),mean = c(m[1],-m[2]),sigma = Gammaneg)$Estimation
  #norm part
  ddnorm1 = dnorm(0,m[1],sqrt(Gamma[1,1]))
  ddnorm2 = dnorm(0,m[2],sqrt(Gamma[2,2]))
  tilmi     = m[2] - Gamma[1,2]*m[1]/Gamma[1,1]
  tilmj     = m[1] - Gamma[1,2]*m[2]/Gamma[2,2]
  tilgammai = Gamma[2,2] - Gamma[1,2]^2/Gamma[1,1]
  tilgammaj = Gamma[1,1] - Gamma[1,2]^2/Gamma[2,2]
  ppnorm1 = pnorm(0,tilmi,sqrt(tilgammai))
  ppnorm2 = pnorm(0,tilmj,sqrt(tilgammaj))

  #for i
  tilSj     = Sigma[-1,-1] - (Sigma[-1,1]/Sigma[1,1])%*%t(Sigma[1,-1])
  cj        = (1 + t(varphi[-1])%*%tilSj%*%varphi[-1])^(-0.5)
  tilvarphi = varphi[1] + (Sigma[-1,1]/Sigma[1,1])%*%varphi[-1]
  #
  tilphi1 = dmvESN1(0,mu[1],Sigma = Sigma[1,1],lambda = sqrt(Sigma[1,1])*cj*tilvarphi,tau = cj*tau)
  tilPhi1 = AcumESN(0,mu[-1] - Sigma[-1,1]/Sigma[1,1]*mu[1],tilSj,c(sqrtm(tilSj)%*%varphi[-1]),c(tau - tilvarphi*mu[1]))

  #for j
  tilSj     = Sigma[-2,-2] - (Sigma[-2,2]/Sigma[2,2])%*%t(Sigma[2,-2])
  cj        = (1 + t(varphi[-2])%*%tilSj%*%varphi[-2])^(-0.5)
  tilvarphi = varphi[2] + (Sigma[-2,2]/Sigma[2,2])%*%varphi[-2]
  #
  tilphi2 = dmvESN1(0,mu[2],Sigma = Sigma[2,2],lambda = sqrt(Sigma[2,2])*cj*tilvarphi,tau = cj*tau)
  tilPhi2 = AcumESN(0,mu[-2] - Sigma[-2,2]/Sigma[2,2]*mu[2],tilSj,c(sqrtm(tilSj)%*%varphi[-2]),c(tau - tilvarphi*mu[2]))

  Eij = (prod(mu) + Sigma[1,2])*(1-2*(PhiESN1+PhiESN2)) +
    (delta[1]*mu[2] + delta[2]*m[1])*(1-2*(Phi1+Phi2)) +
    2*mu[2]*(Sigma[1,1]*tilphi1*(1 - 2*tilPhi1) + Sigma[1,2]*tilphi2*(1 - 2*tilPhi2)) +
    2*delta[2]*(Gamma[1,1]*ddnorm1*(1 - 2*ppnorm1) + Gamma[1,2]*ddnorm2*(1 - 2*ppnorm2)) +
    2*Sigma[2,2]*tilphi2*onlymeanFESN_uni(mu = mu[-2] - Sigma[-2,2]/Sigma[2,2]*mu[2],Sigma = c(tilSj),
                                          lambda = c(sqrtm(tilSj)%*%varphi[-2]),tau = c(tau - tilvarphi*mu[2]))
  return(Eij)
}
