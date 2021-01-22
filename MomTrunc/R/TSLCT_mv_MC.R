# xi = mu
# Omega = Sigma
# lower_p = c(-3,-2,-1)
# upper_p = c(3,4,2)
#
# lower_q = rep(0,length(xi)-length(lower_p))
# upper_q = rep(Inf,length(xi)-length(lower_p))

MCmeanvarTSLCT = function(n = 5000,lower_p,upper_p,xi,Omega,nu=NULL,dist,lower_q,upper_q){
  p = length(c(lower_p))
  q = length(c(lower_q))
  if(length(upper_p) != p)stop("Upper_p dimension does not match lower_p dimension.")
  if(length(upper_q) != q)stop("Upper_p dimension does not match lower_q dimension.")
  if(length(xi) != p+q | ncol(Omega) != p+q | nrow(Omega) != p+q)stop("Xi and Omega with non conformable dimensions. See manual.")
  if(dist != "normal" | dist != "t")stop("The dist values are 'normal' and 't'.")
  if(any(is.na(c(lower_p,upper_p,lower_q,upper_q))))stop("Check limits lower and upper. NA's have been found.")
  res = MCmeanvarTMD(n,lower = c(lower_q,lower_p),upper = c(upper_q,upper_p),mu = xi,Sigma = Omega,nu = nu,dist = dist)
  drop = -(1:q)
  res$mean = as.matrix(res$mean[drop,])
  res$EYY = as.matrix(res$EYY[drop,drop])
  res$varcov = as.matrix(res$varcov[drop,drop])
  return(res)
}


MCmeanvarTSLCT0 = function(n = 5000,lower_p,upper_p,xi,Omega,nu=NULL,dist,lower_q,upper_q){
  q = length(c(lower_q))
  if(dist == "normal"){
    
    res = RcppMCT.lin(n = n,a = c(lower_q,lower_p),b = c(upper_q,upper_p),mu = xi,S = Omega)
  }
  if(dist == "t"){
    if(is.null(nu)){
      stop("Degrees of freedom 'nu' must be provided for the T case.")
    }else{
      if(nu<=0){
        stop("Degrees of freedom 'nu' must be a positive number.")
      }else{
        if(nu >= 300){
          #warning("For degrees of freedom >= 300, Normal case is considered.",immediate. = TRUE)
          res = RcppMCT.lin(n = n,a = c(lower_q,lower_p),b = c(upper_q,upper_p),mu = xi,S = Omega)
        }else{
          res = RcppMCT.lin(n = n,a = c(lower_q,lower_p),b = c(upper_q,upper_p),mu = xi,S = Omega,nu = nu)
        }
      }
    }
  }
  drop = -(1:q)
  res$mean = as.matrix(res$mean[drop,])
  res$EYY = as.matrix(res$EYY[drop,drop])
  res$varcov = as.matrix(res$varcov[drop,drop])
  return(res)
}