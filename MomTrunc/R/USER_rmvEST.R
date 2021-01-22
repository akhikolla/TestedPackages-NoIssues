#######################################################################################
#######################################################################################

rmvEST<-function(n,mu=rep(0,length(lambda)),Sigma=diag(length(lambda)),lambda,tau,nu=NULL){
  if(is.null(nu))stop("Degrees of freedom 'nu' must be provided.")
  #Validating Lambda
  if(nu > 300){return(rmvESN(n,mu,Sigma,lambda,tau))}
  if(is.null(lambda)){
    #not provided by user
    stop("Skewness parameter 'lambda' must be provided (zero vector for the symmetric case).")
  }else{
    #validate input
    if(length(c(lambda)) != length(c(mu)) | !is.numeric(lambda))stop("Lambda must be numeric and have same dimension than mu.")
    if(all(lambda==0)){
      warning("Lambda = 0, symmetric T case is considered.",immediate. = TRUE)
      out = rmvt(n = n,sigma = Sigma,df = nu,delta = mu,type = "shifted")
    }
  }
  if(is.null(tau)){
    #not provided by user
    stop("Extension parameter 'tau' must be provided for the EST case (zero for the Skew-t case).")
  }else{
    #validate input
    if(!is.numeric(tau) | length(tau)>1)stop("Tau must be numeric real number.")
    if(tau == 0){
      return(rmvST(n,mu,Sigma,lambda,nu))
    }
    tautil = tau/sqrt(1+sum(lambda^2))
    if(pt(tautil,nu)< pnorm(-37)){
      #print("normal aproximation")
      Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
      omega_tau = (nu+tautil^2)/(nu+1)
      return(rmvt(n = n,delta = c(mu - tautil*Delta),sigma = omega_tau*(Sigma - Delta%*%t(Delta)),df = nu+1,type = "shifted"))
    }
    return(rEST0(n = n,mu = mu,Sigma = Sigma,lambda = lambda,tau = tau,nu = nu))
  }
}

rEST0<-function(n = 10000,mu=c(0,2),Sigma=diag(2),lambda=c(-1,3),tau=1,nu=4){
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
  #algo = ifelse(pt(tautil,nu)< 10^-3,"gibbs","rejection")
  
  return(RcppTT.GS(n=n,mu = c(mu,0),S = Omega,nu = nu,upper=c(rep(Inf,p),tautil))[,1:p])
}

#meanvarTMD(lower = c(-Inf,-Inf),upper = c(Inf,Inf),mu = c(mu),Sigma = Sigma,lambda = c(lambda),tau = tau,nu = nu,dist = "EST")