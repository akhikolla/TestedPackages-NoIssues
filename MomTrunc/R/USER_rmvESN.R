#######################################################################################
#######################################################################################

rmvESN<-function(n,mu=rep(0,length(lambda)),Sigma=diag(length(lambda)),lambda,tau=0){
  #Validating Lambda
  if(is.null(lambda)){
    #not provided by user
    stop("Skewness parameter 'lambda' must be provided (zero vector for the symmetric case).")
  }else{
    #validate input
    if(length(c(lambda)) != length(c(mu)) | !is.numeric(lambda))stop("Lambda must be numeric and have same dimension than mu.")
    if(all(lambda==0)){
      warning("Lambda = 0, Normal case is considered.",immediate. = TRUE)
      out = rmvnorm(n = n,mean = mu,sigma = Sigma)
    }
  }
  if(is.null(tau)){
    #not provided by user
    stop("Extension parameter 'tau' must be provided for the ESN case (zero for the Skew-normal case).")
  }else{
    #validate input
    if(!is.numeric(tau) | length(tau)>1)stop("Tau must be numeric real number.")
    if(tau == 0){
      return(rmvSN0(n,mu,Sigma,lambda))
    }
    tautil = tau/sqrt(1+sum(lambda^2))
    if(tautil< -37){
      #print("normal aproximation")
      Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
      return(rmvnorm(n = n,mean = c(mu - tautil*Delta),sigma = Sigma - Delta%*%t(Delta)))
    }
    return(rESN0(n = n,mu = mu,Sigma = Sigma,lambda = lambda,tau = tau))
  }
}

rESN0<-function(n = 10000,mu=c(0,2),Sigma=diag(2),lambda=c(-1,3),tau=1){
  mu<-as.matrix(mu)
  Sigma = as.matrix(Sigma)
  lambda<-as.matrix(lambda)
  varphi<-lambda/sqrt(1+sum(lambda^2))
  tautil<-tau/sqrt(1+sum(lambda^2))
  p<-length(mu)
  SS = sqrtm(Sigma)
  Omega1<- cbind(Sigma,-SS%*%varphi)
  Omega2<- cbind(-t(SS%*%varphi),1)
  Omega<- rbind(Omega1,Omega2)
  #algo = ifelse(pnorm(tautil)< 10^-3,"gibbs","rejection")
  return(RcppTT.GS(n=n,mu = c(mu,0),S = Omega,nu = 1000,upper=c(rep(Inf,p),tautil))[,1:p])
}

# gen0 = rtmvnorm(n = n,mean = c(mu,0),sigma = Omega,upper = c(rep(Inf,p),tautil),algorithm = algo)[,1:p]
# colMeans(gen0)
# 
# gen1 = MomTrunc::RcppTT.GS(n=n,mu = c(mu,0),S = Omega,nu = 1000,upper=c(rep(Inf,p),tautil))[,1:p]
# colMeans(gen1)


# rmvSN = function(n,mu,Sigma,lambda){
#   p = length(lambda)
#   T = abs(rnorm(n))
#   Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
#   Gamma = Sigma - Delta%*%t(Delta)
#   gen  = matrix(NA,n,p)
#   for(i in 1:n){
#     gen[i,] = as.matrix(mu) + T[i]*Delta + sqrtm(Gamma)%*%as.matrix(rnorm(p))
#   }
#   return(gen)
# }