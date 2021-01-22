#######################################################################################
#######################################################################################

pmvESN = function(lower = rep(-Inf,length(lambda)),upper=rep(Inf,length(lambda)),mu = rep(0,length(lambda)),Sigma,lambda,tau,log2 = FALSE){
  tautil<-tau/sqrt(1+sum(lambda^2))
  if(tautil< -37){
    #print("normal aproximation")
    Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
    Gamma = Sigma - Delta%*%t(Delta)
    rownames(Gamma) <- colnames(Gamma)
    return(pmvnormt(lower = lower,upper = upper,mean = c(mu - tautil*Delta),sigma = Gamma,uselog2 = log2))
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
    return(pmvnormt(lower = aaum,upper = baum,mean = rep(0,p+1),sigma = Omega,uselog2 = TRUE) - pnorm(tautil,log.p = TRUE)/log(2))
  }else{
    return(pmvnormt(lower = aaum,upper = baum,mean = rep(0,p+1),sigma = Omega)/pnorm(tautil))
  }
}
