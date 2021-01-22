#FULL MOMENTS FOLDED

momentsFMD = function(kappa,mu,Sigma,lambda = NULL,tau = NULL,nu = NULL,dist)
{
  mu = c(mu)
  lambda = c(lambda)
  
  if(!all(c(is.finite(mu)),c(is.finite(Sigma)))){stop("mu and Sigma must contain only finite values.")}
  
  #Validating dims data set
  
  if(dist == "normal" | dist == "SN" | dist == "ESN"){
    if(ncol(as.matrix(kappa)) > 1 | !all(kappa >= 0) | length(c(kappa)) != length(c(mu))) stop("kappa must be a non zero integer vector with same dimensions than mu.")
  }
  
  if(dist == "t" | dist == "ST" | dist =="EST"){
    if(ncol(as.matrix(kappa)) > 1 | !all(kappa >= 0) | !(length(c(kappa)) == length(c(mu)) | length(c(kappa)) == 1)) stop("kappa must be either a non zero integer vector with same dimensions than mu or an scalar being sum(kappa).")
    if(length(c(kappa)) != 1){kappa = sum(kappa)}
  }
  
  if(ncol(as.matrix(mu)) > 1 | !is.numeric(mu)) stop("mu must be numeric and have just one column")
  
  #validate mean an Sigma dimensions
  
  if(ncol(as.matrix(Sigma)) != length(c(mu)))stop("Unconformable dimensions between mu and Sigma")
  if(length(Sigma) == 1){
    if(c(Sigma)<=0)stop("Sigma (sigma^2 for p = 1) must be positive.")
  }else{
    if(!is.pd(Sigma))stop("Sigma must be a square symmetrical real posite definite matrix.")
  }
  
  #validating distributions and nu parameter
  if(dist=="normal"){
    out = RcppKmomentFN(k = kappa,mu = mu,Sigma = Sigma)
  }else
  {
    if(dist == "t"){
      if(is.null(nu)){
        stop("Degrees of freedom 'nu' must be provided for the T case.")
      }else
      {
        if(nu <= sum(kappa)){
          stop("Degrees of freedom 'nu' must be an integer greater than kappa.")
        }else
        {
          if(nu >= 300){
            #warning("For degrees of freedom >= 300, Normal case is considered.",immediate. = TRUE)
            out = RcppKmomentFN(k = kappa,mu = mu,Sigma = Sigma)
          }else
          {
            out = RcppKmomentFT(k = kappa,mu = mu,Sigma = Sigma,nu = nu)
          }
        }
      }
    }else
    {
      if(dist == "ESN" | dist == "SN"){
        #Validating Lambda
        if(is.null(lambda)){
          #not provided by user
          stop("Skewness parameter 'lambda' must be provided for the ESN/SN case.")
        }else
        {
          #validate input
          if(length(c(lambda)) != length(c(mu)) | !is.numeric(lambda))stop("Lambda must be numeric and have same dimension than mu.")
          if(all(lambda==0)){
            warning("Lambda = 0, Normal case is considered.",immediate. = TRUE)
            out = RcppKmomentFN(k = kappa,mu = mu,Sigma = Sigma)
          }
        }
        if(dist=="SN"){
          out = KmomentFESN(k = kappa,mu = mu,Sigma = Sigma,lambda = lambda,tau = 0)
        }else
        {
          if(is.null(tau)){
            #not provided by user
            stop("Extension parameter 'tau' must be provided for the ESN case.")
          }else
          {
            #validate input
            if(!is.numeric(tau) | length(tau)>1)stop("Tau must be numeric real number.")
            out = KmomentFESN(k = kappa,mu = mu,Sigma = Sigma,lambda = lambda,tau = tau)
          }
        }
      }else
      {
        stop("The dist values are 'normal', 't', 'SN' and 'ESN'.")
      }
    }
  }
  # cat('\n')
  # call <- match.call()
  # cat("Call:\n")
  # print(call)
  # cat('\n')
  # print(out)
  # cat('\n')
  return(out)
}
