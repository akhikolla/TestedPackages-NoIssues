dmvEST = function(x,mu=rep(0,length(lambda)),Sigma=diag(length(lambda)),lambda,tau=0,nu){
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
      warning("Lambda = 0, t case is considered.",immediate. = TRUE)
      return(dmvt(x = x,delta = c(mu),sigma = as.matrix(Sigma),df = nu,type = "shifted",log = FALSE))
    }
  }
  if(is.null(tau)){
    #not provided by user
    stop("Extension parameter 'tau' must be provided for the EST case (zero for the Skew-t case).")
  }else{
    #validate input
    if(!is.numeric(tau) | length(tau)>1){
      stop("Tau must be numeric real number.")
    }
    tautil = tau/sqrt(1+sum(lambda^2))
    if(pt(tautil,nu)< pnorm(-37)){
      Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
      Gamma = Sigma - Delta%*%t(Delta)
      rownames(Gamma) <- colnames(Gamma)
      omega_tau = (nu+tautil^2)/(nu+1)
      return(dmvt(x = x,delta = c(mu - tautil*Delta),sigma = omega_tau*Gamma,df = nu+1,type = "shifted",log = FALSE))
    }
    return(dmvEST0(x,mu,Sigma,lambda,tau,nu))
  }
}

########################################################################################################################
########################################################################################################################

dmvEST0 <- function(y, mu, Sigma, lambda,tau,nu){
  #y: deve ser uma matrix onde cada linha tem um vetor de dados multivariados de dimens?o ncol(y) = p. nrow(y) = tamanho da amostra
  #mu, lambda: devem ser do tipo vetor de mesma dimens?o igual a ncol(y) = p
  #Sigma: Matrix p x p
  if(length(c(mu)) == 1){return(dmvEST1(c(y),mu,Sigma,lambda,tau,nu))}
  if(!is.matrix(y)){
    y = matrix(y,ncol = nrow(Sigma),byrow = TRUE)
  }
  n <- nrow(y)
  p <- ncol(y)
  tautil<-tau/sqrt(1+sum(lambda^2))
  
  nu2y = (nu + p)/(nu + mahalanobis(y,center = mu,cov = Sigma))
  
  dens <- dmvt(x = y,delta = c(mu),sigma = Sigma,df = nu,type = "shifted",log = FALSE)*exp(
    pt(sqrt(nu2y)*(
         apply(matrix(rep(t(lambda)%*%solve(sqrtm(Sigma)),n), n, p, byrow = TRUE)*
               (y - matrix(rep(mu, n), n, p, byrow = TRUE)), 1,sum) + tau),
       df = nu+p,log.p = TRUE) - pt(tautil,nu,log.p = TRUE))
  return(dens)
}

########################################################################################################################
########################################################################################################################

dmvEST1 <- function(y, mu=0, Sigma=1, lambda,tau,nu){
  #y: deve ser uma matrix onde cada linha tem um vetor de dados multivariados de dimens?o ncol(y) = p. nrow(y) = tamanho da amostra
  #mu, lambda: devem ser do tipo vetor de mesma dimens?o igual a ncol(y) = p
  #Sigma: Matrix p x p
  n <- length(c(y))
  p <- 1
  s = sqrt(Sigma)
  tautil = tau/sqrt(1+sum(lambda^2))
  if(pt(tautil,nu)< pnorm(-35)){
    #print("normal aproximation")
    Gamma  = Sigma/(1+lambda^2)
    mub    = lambda*tau*Gamma/s
    omega_tau = (nu+tautil^2)/(nu+1)
    return(dent(y,mu-mub,omega_tau*Gamma,nu+1))
  }
  
  nu2y = (nu + 1)/(nu + mahalanobis(y,center = mu,cov = Sigma))
  
  dens <- dent(x = c(y),mu = c(mu),sigma2 = Sigma,nu = nu)*exp(
    pt(
      sqrt(nu2y)*(
      apply(matrix(rep(t(lambda)%*%solve(sqrtm(Sigma)),n), n, p, byrow = TRUE)*
               (y - matrix(rep(mu, n), n, p, byrow = TRUE)), 1,sum) + tau),
      df = nu+1,log.p = TRUE) - pt(tautil,nu,log.p = TRUE))
  return(dens)
}