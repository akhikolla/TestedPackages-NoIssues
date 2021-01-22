
#location-scale student-T distribution pdf
dent<-function(x,mu,sigma2,nu,log = FALSE){
  z<-(x-mu)/sqrt(sigma2)
  if(log){
    return(log(1/sqrt(sigma2)) + dt(x = z,df = nu,log = TRUE))
  }else{
  return(1/sqrt(sigma2)*dt(x = z,df = nu))
    }
}

#########################################################################

#location-scale student-T distribution pcf
pent<-function(x,mu,sigma2,nu,log = FALSE){
  return(pt((x-mu)/sqrt(sigma2),nu,log.p = log))
}

#########################################################################
#
#######################################################################################
#######################################################################################


pent2 = function(lower = -Inf,upper = Inf,mu = 0,sigma2 = 1,nu,...){
  if(all(lower == -Inf)){
    return(pt((upper-mu)/sqrt(sigma2),df = nu,...))
  }else{
    if(all(upper == Inf)){
      return(pt((lower-mu)/sqrt(sigma2),df = nu,lower.tail = FALSE,...))
    }else{
      return(pt((upper-mu)/sqrt(sigma2),nu) - pt((lower-mu)/sqrt(sigma2),nu))
    }
  }
}

#######################################################################################
#######################################################################################

invmillsT = function(x,mu=0,sigma2=1,nu){
  return(exp(dent(x,mu,sigma2,nu,log = TRUE) - pent2(upper = x,mu = mu,sigma2 = sigma2,nu = nu,log.p = TRUE)))
}
