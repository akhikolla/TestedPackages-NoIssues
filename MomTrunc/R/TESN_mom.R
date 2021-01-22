KmomentESN = function(k,a,b,mu,Sigma,lambda,tau)
{
  #tau goes to infinite
  tautil<-tau/sqrt(1+sum(lambda^2))
  if(tautil< -36){
    #print("normal aproximation")
    Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
    return(KmomentN(k = k,a = a, b = b,mu = c(mu - tautil*Delta),Sigma = Sigma - Delta%*%t(Delta)))
  }
  p = length(mu)
  SS   = sqrtm(Sigma)
  tautil = tau/sqrt(1+sum(lambda^2))
  varpsi = lambda/sqrt(1+sum(lambda^2))
  Omega  = cbind(rbind(Sigma,-t(varpsi)%*%SS),rbind(-SS%*%varpsi,1))
  rownames(Omega) <- colnames(Omega)
  return(KmomentN(c(k,0),a = c(a,-10^7),b = c(b,tautil),mu = c(mu,0),Sigma = Omega)[,-(p+1)])
}
