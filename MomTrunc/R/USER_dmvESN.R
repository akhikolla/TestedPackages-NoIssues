dmvESN = function(x,mu=rep(0,length(lambda)),Sigma=diag(length(lambda)),lambda,tau=0){
  #Validating Lambda
  if(is.null(lambda)){
    #not provided by user
    stop("Skewness parameter 'lambda' must be provided (zero vector fot the symmetric case).")
  }else{
    #validate input
    if(length(c(lambda)) != length(c(mu)) | !is.numeric(lambda)){
      stop("Lambda must be numeric and have same dimension than mu.")
    }
    if(all(lambda==0)){
      warning("Lambda = 0, Normal case is considered.",immediate. = TRUE)
      return(dmvnorm(x = x,mean = c(mu),sigma = as.matrix(Sigma)))
    }
  }
  if(is.null(tau)){
    #not provided by user
    stop("Extension parameter 'tau' must be provided for the ESN case (zero for the Skew-normal case).")
  }else{
    #validate input
    if(!is.numeric(tau) | length(tau)>1){
      stop("Tau must be numeric real number.")
    }
    tautil = tau/sqrt(1+sum(lambda^2))
    if(tautil< -37){
      #print("normal aproximation")
      Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
      return(dmvnorm(x = x,mean = c(mu - tautil*Delta),sigma = as.matrix(Sigma - Delta%*%t(Delta))))
    }
    return(dmvESN0(x,mu,Sigma,lambda,tau))
  }
}

dmvESN0 <- function(y, mu, Sigma, lambda,tau){
  #y: deve ser uma matrix onde cada linha tem um vetor de dados multivariados de dimens?o ncol(y) = p. nrow(y) = tamanho da amostra
  #mu, lambda: devem ser do tipo vetor de mesma dimens?o igual a ncol(y) = p
  #Sigma: Matrix p x p
  if(length(c(mu)) == 1){return(dmvESN1(c(y),mu,Sigma,lambda,tau))}
  if(!is.matrix(y)){
    y = matrix(y,ncol = nrow(Sigma),byrow = TRUE)
  }
  n <- nrow(y)
  p <- ncol(y)
  tautil<-tau/sqrt(1+sum(lambda^2))
  dens <- dmvnorm(y, c(mu),Sigma)*exp(
    pnorm(apply(matrix(rep(t(lambda)%*%solve(sqrtm(Sigma)),n), n, p, byrow = TRUE)*
                  (y - matrix(rep(mu, n), n, p, byrow = TRUE)), 1,sum)+tau,log.p = TRUE) - pnorm(tautil,log.p = TRUE))
  return(dens)
}

#######################################################################################
#######################################################################################

dmvESN1 <- function(y, mu=0, Sigma=1, lambda,tau){
  #y: deve ser uma matrix onde cada linha tem um vetor de dados multivariados de dimens?o ncol(y) = p. nrow(y) = tamanho da amostra
  #mu, lambda: devem ser do tipo vetor de mesma dimens?o igual a ncol(y) = p
  #Sigma: Matrix p x p
  n <- length(c(y))
  p <- 1
  s = sqrt(Sigma)
  tautil = tau/sqrt(1+sum(lambda^2))
  if(tautil< -35){
    #print("normal aproximation")
    Gamma  = Sigma/(1+lambda^2)
    mub    = lambda*tau*Gamma/s
    return(dnorm(y,mu-mub,sqrt(Gamma)))
  }
  dens <- dnorm(x = c(y),mean = c(mu),sd = sqrt(Sigma))*exp(
      pnorm(apply(matrix(rep(t(lambda)%*%solve(sqrtm(Sigma)),n), n, p, byrow = TRUE)*
                  (y - matrix(rep(mu, n), n, p, byrow = TRUE)), 1,sum)+tau,log.p = TRUE) - pnorm(tautil,log.p = TRUE))
  return(dens)
}
