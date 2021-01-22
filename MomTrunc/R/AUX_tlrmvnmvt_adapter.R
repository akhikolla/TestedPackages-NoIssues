pmvn.genz = function(lower = rep(-Inf,ncol(sigma)),upper = rep(Inf,ncol(sigma)),mean = rep(0,ncol(sigma)),sigma,uselog2 = FALSE){
  
  return(list(Estimation = tlrmvnmvt::pmvn(lower = lower,upper = upper,mean = mean,sigma = sigma,uselog2)[1]))
  
}

# pmvn.genz(lower = lower,upper = upper,mean = mu,sigma = Sigma,uselog2 = FALSE)$Estimation


pmvt.genz = function(lower = rep(-Inf,ncol(sigma)),upper = rep(Inf,ncol(sigma)),mean = rep(0,ncol(sigma)),sigma,nu,uselog2 = FALSE,N = NULL){
  
  if(is.null(N)){
    return(list(Estimation = tlrmvnmvt::pmvt(lower = lower,upper = upper,delta = mean,df = nu,sigma = sigma,uselog2,type = "shifted")[1]))
  }else{
    return(list(Estimation = tlrmvnmvt::pmvt(lower = lower,upper = upper,delta = mean,df = nu,sigma = sigma,uselog2,algorithm = GenzBretz(N = N),type = "shifted")[1]))
  }
}

# pmvt.genz(lower = a,upper = b,mean = mu,nu = nu,sigma = Sigma,uselog2 = TRUE,N = 799)$Estimation
# 
# pmvt.genz(lower = a,upper = b,mean = mu,nu = nu,sigma = Sigma,uselog2 = TRUE)
# 
# pmvt.genz(lower = a,upper = b,mean = mu,nu = nu,sigma = Sigma)
