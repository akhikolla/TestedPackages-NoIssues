onlymeanTSLCT0 = function(lower_p,upper_p,xi,Omega,nu=NULL,dist,lower_q,upper_q){
  q = length(c(lower_q))
  if(dist == "normal"){
    res = onlymeanN(lower = c(lower_q,lower_p),upper = c(upper_q,upper_p),mu = xi,Sigma = Omega)
  }
  if(dist == "t"){
    if(is.null(nu)){
      stop("Degrees of freedom 'nu' must be provided for the T case.")
    }else{
      if(nu%%1!=0){
        stop("Degrees of freedom 'nu' must be an integer greater than 2.")
      }else{
        if(nu <= 2){stop("Sorry, we can only compute the first moment for degrees of freedom larger than 2.")
        }else{
          if(nu >= 200){
            warning("For degrees of freedom >= 200, Normal case is considered.",immediate. = TRUE)
            res = onlymeanN(lower = c(lower_q,lower_p),upper = c(upper_q,upper_p),mu = xi,Sigma = Omega)
          }else{
            if(nu < 4){
              warning("Sorry, we can only compute the second moment when the degrees of freedom is larger than 3.",immediate. = TRUE)
              res = onlymeanTall(lower = c(lower_q,lower_p),upper = c(upper_q,upper_p),mu = xi,Sigma = Omega,nu = nu)
            }else{
              res = onlymeanTall(lower = c(lower_q,lower_p),upper = c(upper_q,upper_p),mu = xi,Sigma = Omega,nu = nu)
            }
          }
        }
      }
    }
  }
  drop = -(1:q)
  res$mean = as.matrix(res$mean[drop,])
  return(res)
}
