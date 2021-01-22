optmeanvarTMD = function(lower = rep(-Inf,length(mu)),upper = rep(Inf,length(mu)),mu,Sigma,lambda = NULL, tau = NULL,Gamma = NULL,nu = NULL,dist,n = 10000){
  
  p = length(mu)
  
  bool1 = is.infinite(lower)  
  
  bool2 = is.infinite(upper)  
  
  infinite.dims = sum(bool1*bool2)
  
  dim = length(c(mu)) + ifelse(is.null(lambda),0,ncol(as.matrix(lambda))) - infinite.dims
  
  if(infinite.dims < p){
    
    if(all(is.infinite(lower[!(bool1 & bool2)])) | all(is.infinite(upper[!(bool1 & bool2)]))){
      
      if(dist %in% c("normal","SN","ESN","SUN")){
        
        dim = dim - 2}
        
      # }else{
      #   
      #   dim = dim - 1
      #   
      # }
      
    }
    
  }else{
    
    return(MomTrunc::meanvarTMD(lower,upper,mu,Sigma,lambda,tau,Gamma,nu,dist))
    
  }
  
  if(dim < 3){
    
    return(MomTrunc::meanvarTMD(lower,upper,mu,Sigma,lambda,tau,Gamma,nu,dist))
    
  }else{
    
    return(MomTrunc::MCmeanvarTMD(lower,upper,mu,Sigma,lambda,tau,Gamma,nu,dist,n))
    
  }
  
}

