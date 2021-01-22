#MEAN AND VARIANCE

meanvarFMD = function(mu,Sigma,lambda = NULL,tau = NULL,nu = NULL,dist)
{
  mu = c(mu)
  lambda = c(lambda)

  if(!all(c(is.finite(mu)),c(is.finite(Sigma)))){stop("mu and Sigma must contain only finite values.")}

  #Validating dims data set
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
    out = RcppmeanvarFN(mu = mu,S = Sigma)
  }else{
    if(dist == "t"){
      if(is.null(nu)){
        stop("Degrees of freedom 'nu' must be provided for the T case.")
      }else{
        if(nu%%1!=0){
          stop("Degrees of freedom 'nu' must be an integer greater than 1.")
        }else{
          if(nu <= 1){stop("Sorry, moments only exists for degrees of freedom larger than 1.")
          }else{
            if(nu >= 300){
              #warning("For degrees of freedom >= 300, Normal case is considered.",immediate. = TRUE)
              out = RcppmeanvarFN(mu = mu,S = Sigma)
            }else{
              if(nu == 2){
                warning("Second moment does not exist for degrees of freedom less than 3.",immediate. = TRUE)
                out = RcppmeanvarFT(mu = mu,S = Sigma,nu = nu)
              }else{
                out = RcppmeanvarFT(mu = mu,S = Sigma,nu = nu)
              }
            }
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
            out = RcppmeanvarFN(mu = mu,S = Sigma)
          }
        }
        if(dist=="SN"){
          out = meanvarFESNopt(mu = mu,Sigma = Sigma,lambda = lambda,tau = 0)
        }else{
          if(is.null(tau)){
            #not provided by user
            stop("Extension parameter 'tau' must be provided for the ESN case.")
          }else{
            #validate input
            if(!is.numeric(tau) | length(tau)>1)stop("Tau must be numeric real number.")
            out = meanvarFESNopt(mu = mu,Sigma = Sigma,lambda = lambda,tau = tau)
          }
        }
      }else{
        stop("The dist values are 'normal', 't', 'SN' and 'ESN'.")
      }
    }
  }
  # cat('\n')
  # call <- match.call()
  # cat("Call:\n")
  # print(call)
  # cat('\n')
  # cat("Mean:\n")
  # print(out$mean)
  # if(dist == "normal" | (dist == "t" & nu >= 4)){
  #   cat('\n')
  #   if(length(mu)==1){cat("Variance:\n")}else{cat("Varcov matrix:\n")}
  #   print(out$varcov)
  #   cat('\n')
  # }
  return(out)
}

#TESTING
# a = c(-0.8,-0.7,-0.6,-0.5)
# b = c(0.5,0.6,0.7,0.8)
# mu = c(0.1,0.2,0.3,0.4)
# S = matrix(data = c(1,0.2,0.3,0.1,0.2,1,0.4,-0.1,0.3,0.4,1,0.2,0.1,-0.1,0.2,1),nrow = length(mu),ncol = length(mu),byrow = TRUE)

# resN = meanvarFMD(mu,S,dist = "normal")
# resT = meanvarFMD(mu,S,dist = "t",nu=4)
