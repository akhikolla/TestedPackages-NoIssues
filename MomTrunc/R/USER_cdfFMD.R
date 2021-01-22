cdfFMD = function(x,mu,Sigma,lambda = NULL,tau = NULL,dist,nu = NULL)
{
  #Validating dims data set
  if(ncol(as.matrix(mu)) > 1 | !is.numeric(mu)) stop("y must be numeric and have just one column")
  if(ncol(as.matrix(x)) > 1 | !all(x >= 0) | length(c(x)) != length(c(mu))) stop("x must be numeric, positive, with same dimensions than mu.")

  #validate mean an Sigma dimensions

  if(ncol(as.matrix(Sigma)) != length(c(mu)))stop("Unconformable dimensions between mu and Sigma")
  if(length(Sigma) == 1){
    if(c(Sigma)<=0)stop("Sigma (sigma^2 for p = 1) must be positive.")
  }else{
    if(!is.pd(Sigma))stop("Sigma must be a square symmetrical real posite definite matrix.")
  }

  #validating distributions and nu parameter

  if(dist=="normal"){
    out = cdfFN(x = x,mu = mu,Sigma = Sigma)
  }else{
    if(dist == "t"){
      if(is.null(nu)){
        stop("Degrees of freedom 'nu' must be provided for the T case.")
      }else{
        if(nu%%1!=0){
          stop("Degrees of freedom 'nu' must be an integer greater than 2.")
        }else{
          if(nu >= 200){
            warning("For degrees of freedom >= 200, Normal case is considered.",immediate. = TRUE)
            out = cdfFN(x = x,mu = mu,Sigma = Sigma)
          }else{
            out = cdfFT(x = x,mu = mu,Sigma = Sigma,nu = nu)
          }
        }
      }
    }else{
      if(dist == "ESN" | dist == "SN"){
        #Validating Lambda
        if(is.null(lambda)){
          #not provided by user
          stop("Skewness parameter 'lambda' must be provided for the ESN/SN case.")
        }else{
          #validate input
          if(length(c(lambda)) != length(c(mu)) | !is.numeric(lambda))stop("Lambda must be numeric and have same dimension than mu.")
          if(all(lambda==0)){
            warning("Lambda = 0, Normal case is considered.",immediate. = TRUE)
            out = cdfFN(x = x,mu = mu,Sigma = Sigma)
          }
        }
        if(dist=="SN"){
          out = cdfFESN(x,mu,Sigma,lambda,0)
        }else{
          if(is.null(tau)){
            #not provided by user
            stop("Extension parameter 'tau' must be provided for the ESN case.")
          }else{
            #validate input
            if(!is.numeric(tau) | length(tau)>1)stop("Tau must be numeric real number.")
            out = cdfFESN(x,mu,Sigma,lambda,tau)
          }
        }
      }else{
        stop("The dist values are 'normal', 't', 'SN' and 'ESN'.")
      }
    }
  }
  return(out)
}
