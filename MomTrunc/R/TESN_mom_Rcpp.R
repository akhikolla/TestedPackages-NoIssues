RcppKmomentESN = function(k,a,b,mu,Sigma,lambda,tau)
{
  p = length(mu)
  #tau goes to infinite
  tautil<-tau/sqrt(1+sum(lambda^2))
  if(tautil< -36){
    #print("normal aproximation")
    Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
    
    return(RcppKmomentN(k = k,a = a, b = b,mu = mu - tautil*Delta,Sigma = Sigma - Delta%*%t(Delta)))
  }
  SS   = sqrtm(Sigma)
  varpsi = lambda/sqrt(1+sum(lambda^2))
  Omega  = cbind(rbind(Sigma,-t(varpsi)%*%SS),rbind(-SS%*%varpsi,1))
  rownames(Omega) <- colnames(Omega)
  return(RcppKmomentN(c(k,0),a = c(a,-10^7),b = c(b,tautil),mu = c(mu,0),Sigma = Omega)[,-(p+1)])
}


# RcppKmomentESN = function(k,a,b,mu,Sigma,lambda,tau)
# {
#   k = as.matrix(k)
#   a = as.matrix(a)
#   b = as.matrix(b)
#   mu = as.matrix(mu)
#   Sigma = as.matrix(Sigma)
#   
#   #tau goes to infinite
#   tautil<-tau/sqrt(1+sum(lambda^2))
#   if(tautil< -36){
#     #print("normal aproximation")
#     Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
#     
#     return(RcppKmomentN(k = k,a = a, b = b,mu = mu - tautil*Delta,Sigma = Sigma - Delta%*%t(Delta)))
#   }
#   p = length(mu)
#   SS   = sqrtm(Sigma)
#   tautil = tau/sqrt(1+sum(lambda^2))
#   varpsi = lambda/sqrt(1+sum(lambda^2))
#   Omega  = cbind(rbind(Sigma,-t(varpsi)%*%SS),rbind(-SS%*%varpsi,1))
#   rownames(Omega) <- colnames(Omega)
#   return(RcppKmomentN(rbind(k,0),a = rbind(a,-10^7),b = rbind(b,tautil),mu = rbind(mu,0),Sigma = Omega)[,-(p+1)])
# }
