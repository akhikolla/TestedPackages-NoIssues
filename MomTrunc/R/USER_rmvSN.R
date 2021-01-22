#######################################################################################
#######################################################################################

rmvSN<-function(n,mu=rep(0,length(lambda)),Sigma=diag(length(lambda)),lambda){
  #Validating Lambda
  if(is.null(lambda)){
    #not provided by user
    stop("Skewness parameter 'lambda' must be provided (zero vector fot the symmetric case).")
  }else{
    #validate input
    if(length(c(lambda)) != length(c(mu)) | !is.numeric(lambda))stop("Lambda must be numeric and have same dimension than mu.")
    if(all(lambda==0)){
      warning("Lambda = 0, Normal case is considered.",immediate. = TRUE)
      out = rmvnorm(n = n,mean = mu,sigma = Sigma)
    }else{
      return(rmvSN0(n,mu,Sigma,lambda))
    }
  }
}


rmvSN0 = function(n,mu,Sigma,lambda){
  p = length(lambda)
  T = abs(rnorm(n))
  Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
  Gamma = Sigma - Delta%*%t(Delta)
  gen  = matrix(NA,n,p)
  for(i in 1:n){
    gen[i,] = as.matrix(mu) + T[i]*Delta + sqrtm(Gamma)%*%as.matrix(rnorm(p))
  }
  return(gen)
}
