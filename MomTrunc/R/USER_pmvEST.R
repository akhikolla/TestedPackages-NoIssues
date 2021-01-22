#######################################################################################
#######################################################################################
# lower=c(-10,-10)
# upper = c(10,10)
# Sigma=matrix(c(1, -0.5, -0.5, 1), 2, 2)
# lambda = c(1/2,-1/2)
# tau = 1
# nu = 4



pmvEST = function(lower = rep(-Inf,length(lambda)),upper=rep(Inf,length(lambda)),mu = rep(0,length(lambda)),Sigma,lambda,tau,nu,log2 = FALSE){
  tautil<-tau/sqrt(1+sum(lambda^2))
  if(pt(tautil,nu)< pnorm(-37)){
    #print("normal aproximation")
    Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
    Gamma = Sigma - Delta%*%t(Delta)
    rownames(Gamma) <- colnames(Gamma)
    omega_tau = (nu+tautil^2)/(nu+1)
    return(pmvnormt(lower = lower,upper = upper,mean = c(mu - tautil*Delta),sigma = omega_tau*Gamma,nu = nu+1,uselog2 = log2))
  }
  aaum = c(lower-mu,-Inf)
  baum = c(upper-mu,tautil)
  mu<-as.matrix(mu)
  lambda<-as.matrix(lambda)
  varphi<-lambda/sqrt(1+sum(lambda^2))
  p<-length(mu)
  SS = sqrtm(Sigma)
  Omega1<- cbind(Sigma,-SS%*%varphi)
  Omega2<- cbind(-t(SS%*%varphi),1)
  Omega<- rbind(Omega1,Omega2)
  rownames(Omega) <- colnames(Omega)

  if(log2 == TRUE){
    return(min(pmvnormt(lower = aaum,upper = baum,sigma = Omega,nu = nu+1,uselog2 = TRUE) - pt(tautil,nu,log.p = TRUE)/log(2),0))
  }else{
    return(min(pmvnormt(lower = aaum,upper = baum,sigma = Omega,nu = nu+1)/pt(tautil,nu),1))
  }
}