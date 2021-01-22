#This function uses the best (in terms of time/precision) function from the packages above depending of the dimension of the vector, 
#degrees of freedom and even if one method collapses


pmvnormt = function(lower = rep(-Inf,ncol(sigma)),upper = rep(Inf,ncol(sigma)),mean = rep(0,ncol(sigma)),sigma,nu = NULL,uselog2 = FALSE){
  
  mean = c(mean)
  sigma = as.matrix(sigma)
  p = ncol(sigma)
  
  lower = c(lower)
  upper = c(upper)
  
  if(is.null(nu)){
    #normal case
    if(p < 10){
      
      prob = mvtnorm::pmvnorm(lower = lower,upper = upper,mean = mean,sigma = sigma)[1]
      
      if(prob < 0){
        
        return(pmvn.genz(lower = lower,upper = upper,mean = mean,sigma = sigma,uselog2 = uselog2)[[1]])
        
      }else{
        
        return(ifelse(uselog2,log2(prob),prob))
      }
    }else{
      
      return(pmvn.genz(lower = lower,upper = upper,mean = mean,sigma = sigma,uselog2 = uselog2)[[1]])
      
    }
  }else{
    #student t case
    if(p < 10 & nu%%1 == 0){
      lower = lower - mean
      upper = upper - mean
      
      prob = mvtnorm::pmvt(lower = lower,upper = upper,sigma = sigma,df = nu)[1]
      
      if(prob < 0){
        
        return(pmvt.genz(lower = lower,upper = upper,sigma = sigma,nu = nu,uselog2 = uselog2)[[1]])
        
      }else{
        
        return(ifelse(uselog2,log2(prob),prob))
      }
      
    }else{
      return(pmvt.genz(lower = lower,upper = upper,mean = mean,sigma = sigma,nu = nu,uselog2 = uselog2)[[1]])
    }
  }
}
