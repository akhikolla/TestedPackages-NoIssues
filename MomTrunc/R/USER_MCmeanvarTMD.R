#MEAN AND VARIANCE

MCmeanvarTMD = function(lower = rep(-Inf,length(mu)),upper = rep(Inf,length(mu)),mu,Sigma,lambda = NULL, tau = NULL,Gamma = NULL,nu = NULL,dist,n = 10000)
{
  mu = c(mu)
  tau = c(tau)
  
  if(dist == "SN"){
    return(MCmeanvarTMD(n = n,lower = lower,upper = upper,mu = mu,Sigma = Sigma,lambda = lambda,Gamma = 1,tau = 0,dist = "SUN"))
  }
  
  if(dist == "ESN"){
    return(MCmeanvarTMD(n = n,lower = lower,upper = upper,mu = mu,Sigma = Sigma,lambda = lambda,Gamma = 1,tau = tau,dist = "SUN"))
  }
  
  if(dist == "ST"){
    return(MCmeanvarTMD(n = n,lower = lower,upper = upper,mu = mu,Sigma = Sigma,lambda = lambda,Gamma = 1,tau = 0,nu = nu,dist = "SUT"))
  }
  
  if(dist == "EST"){
    return(MCmeanvarTMD(n = n,lower = lower,upper = upper,mu = mu,Sigma = Sigma,lambda = lambda,Gamma = 1,tau = tau,nu = nu,dist = "SUT"))
  }
  
  #SUN and SUT distributions
  if(dist == "SUN" | dist == "SUT"){
    p = length(mu)
    q = length(tau)
    
    if(!all(c(is.finite(mu)),c(is.finite(Sigma)))){stop("mu and Sigma must contain only finite values.")}
    
    if(any(is.na(c(mu,Sigma,lambda,tau,Gamma))))stop("Check parameters mu, Sigma, lambda, tau and Gamma. NA's have been found.")
    
    if(is.null(lambda) | is.null(tau) | is.null(Gamma)) stop("Lambda, tau and Gamma parameters must be provided.")
    
    lambda = as.matrix(lambda)
    Gamma = as.matrix(Gamma)
    
    if(nrow(Sigma) != p | ncol(Sigma) != p | nrow(lambda) != p | ncol(lambda) != q | nrow(Gamma) != q | ncol(Gamma) != q){
      stop("Nonconforming parameters dimensions. See manual.")
    }
    
    if(!(all(eigen(Gamma)$values >= 0) & all(diag(Gamma) == 1))){stop("Gamma must be a valid correlation matrix.")}
    
    xi = c(tau,mu)
    Delta = sqrtm(Sigma)%*%lambda
    Omega1 = cbind(Gamma + t(lambda)%*%lambda,t(Delta))
    Omega2 = cbind(Delta,Sigma)
    Omega  = rbind(Omega1,Omega2)
    if(dist == "SUN"){
      out = MCmeanvarTSLCT0(n = n,lower_p = lower,upper_p = upper,xi = xi,Omega = Omega,nu = nu,dist = "normal",lower_q = rep(0,q),upper_q = rep(Inf,q))
    }
    if(dist == "SUT"){
      out = MCmeanvarTSLCT0(n = n,lower_p = lower,upper_p = upper,xi = xi,Omega = Omega,nu = nu,dist = "t",lower_q = rep(0,q),upper_q = rep(Inf,q))
    }
    return(out)
  }
  
  lambda = c(lambda)
  
  if(!all(c(is.finite(mu)),c(is.finite(Sigma)))){stop("mu and Sigma must contain only finite values.")}
  
  #Validating dims data set
  if(ncol(as.matrix(mu)) > 1 | !is.numeric(mu)) stop("mu must be numeric and have just one column")
  
  #validate mean an Sigma dimensions
  if(ncol(as.matrix(Sigma)) != length(c(mu)))stop("Unconformable dimensions between mu and Sigma")
  if(length(Sigma) == 1){
    if(c(Sigma)<=0)stop("Sigma (sigma^2 for p = 1) must be positive.")
  }else{
    if(!is.pd(Sigma))stop("Sigma must be a square symmetrical real positive definite matrix.")
  }
  if(all(is.null(lower))){
    lower = rep(-Inf,length(mu))
  }else{
    if(length(c(lower)) != length(c(mu)) | !is.numeric(lower))stop("Lower bound must be numeric and have same dimension than mu.")
  }
  if(all(is.null(upper))){
    upper = rep(Inf,length(mu))
  }else{
    if(length(c(upper)) != length(c(mu)) | !is.numeric(upper))stop("Upper bound must be numeric and have same dimension than mu.")
  }
  if(all(lower < upper) == FALSE)stop("Lower bound must be lower than or equal to upper bound.")
  
  #validating distributions and nu parameter
  if(dist=="normal"){
    out = RcppMCT.lin(n = n,a = lower,b = upper,mu = mu,S = Sigma)
  }else{
    if(dist == "t"){
      if(is.null(nu)){
        stop("Degrees of freedom 'nu' must be provided for the T case.")
      }else{
        if(nu<=0){
          stop("Degrees of freedom 'nu' must be a positive number.")
        }else{
          if(nu >= 300){
            #warning("For degrees of freedom >= 300, Normal case is considered.",immediate. = TRUE)
            out = RcppMCT.lin(n = n,a = lower,b = upper,mu = mu,S = Sigma)
          }else{
            out = RcppMCT.lin(n = n,a = lower,b = upper,mu = mu,S = Sigma,nu = nu)
          }
        }
      }
    }else{
        stop("The dist values are 'normal', 't', 'SN', 'ST', 'ESN', 'EST', 'SUN' and 'SUT'.")
    }
  }
  return(out)
}
