meanvarTall = function(lower=rep(-Inf,length(mu)),upper=rep(Inf,length(mu)),mu,Sigma,nu,omega = FALSE){
  
  p = length(mu)
  
  if(p == 1){
    
    if(nu >= 3){
      
      return(meanvarT16(a = lower,b = upper,mu = mu,Sigma = Sigma,nu=nu,omega))
      
    }else{
      
      if(omega){

          return(RcppMCT.lin(n = 5000,a = lower,b = upper,mu = mu,S = as.matrix(Sigma),nu = nu,omega = omega))

        
      }else{
          
          return(dtmvtmuvar(a = lower,b = upper,mu = mu,S = Sigma,nu=nu))
        
      }
      
    }
    
  }
  
  
  #Multivariate
  ###########################################################################################################################
  
  bool1 = is.infinite(lower)
  bool2 = is.infinite(upper)
  
  if(sum(bool1*bool2) > 0){ #Does exist infinite pairs?
    
    if(sum(bool1*bool2) == p){ #All infinites?
      
      if(nu > 2){
        
        varcov = nu/(nu-2)*Sigma
        EYY = varcov + mu%*%t(mu)
        
        return(list(mean = mu,EYY = EYY,varcov = varcov))
        
      }else{
        
        return(list(mean = mu,EYY = matrix(NaN,p,p),varcov = matrix(NaN,p,p)))
        
      }
      
    }else{
      
      return(withinfsT(lower,upper,mu,Sigma,nu))
      
    }
  }
  
  
  if(nu > 4){
    
    #Lin algorithm: TTmoment package
    
    if(sum(bool1) + sum(bool2) == 0){ #All limits are finite?
      
      return(meanvarT.Lin.IC(a = lower,b = upper,mu = mu,S = Sigma,nu=nu,omega))
      
    }else{
      
      if(sum(bool1) == p){ #All lower limits are infinite?
        
        
        if(sum(bool2) == p){ #All lower and upper limits are infinite?
          
          #NO TRUNCATION
          
          varcov = nu/(nu-2)*Sigma
          return(list(mean = mu,EYY = varcov + mu%*%t(mu),varcov = varcov))
          
        }else{
          
          #LEFT CENSORING
          
          return(meanvarT.Lin.RC(b = upper,mu = mu,S = Sigma,nu=nu,omega))
        }
        
      }else{
        
        if(sum(bool2) == p){ #All upper limits are infinite?
          
          #RIGHT TRUNCATION
          
          return(meanvarT.Lin.LC(a = lower,mu = mu,S = Sigma,nu=nu,omega))
          
        }else{
          
          #All kind of censoring
          
          return(meanvarT.Lin.LRIC(a = lower,b = upper,mu = mu,S = Sigma,nu=nu,omega))
          
        }
        
      }
      
    }
    
  }
  
  ###########################################################################################################################
  
  if(nu >= 3){
    
    #Galarza et.al. algorithm: T paper
    
    if(sum(bool1) + sum(bool2) == 0){ #All limits are finite?
      
      return(meanvarT16_finite(a = lower,b = upper,mu = mu,Sigma = Sigma,nu=nu,omega))
      
    }else{
      
      if(sum(bool1) == p){ #All lower limits are infinite?
        
        
        if(sum(bool2) == p){ #All lower and upper limits are infinite?
          
          #NO TRUNCATION
          
          varcov = nu/(nu-2)*Sigma
          return(list(mean = mu,EYY = varcov + mu%*%t(mu),varcov = varcov))
          
        }else{
          
          #LEFT CENSORING
          
          return(meanvarT16_upper(b = upper,mu = mu,Sigma = Sigma,nu=nu,omega))
        }
        
      }else{
        
        if(sum(bool2) == p){ #All upper limits are infinite?
          
          #RIGHT TRUNCATION
          
          return(meanvarT16_lower(a = lower,mu = mu,Sigma = Sigma,nu=nu,omega))
          
        }else{
          
          #All kind of censoring
          
          return(meanvarT16(a = lower,b = upper,mu = mu,Sigma = Sigma,nu=nu,omega))
          
        }
        
      }
      
    }
    
  }
  
  ###########################################################################################################################
  
  #nu < 3
  
  #Kan method
  
  if(omega){
    
    #given that we need omega estimation, we must run MC for these cases
      
      return(RcppMCT.lin(n = 5000,a = lower,b = upper,mu = mu,S = as.matrix(Sigma),nu = nu,omega = omega))
    
  }else{
    
    if(sum(bool1) + sum(bool2) == 0){ #All limits are finite?
      
      return(ftmvtmuvar(a = lower,b = upper,mu = mu,S = Sigma,nu=nu))
      
    }else{
      
      if(sum(bool1) == p){ #All lower limits are infinite?
        
        
        if(sum(bool2) == p){ #All lower and upper limits are infinite?
          
          #NO TRUNCATION
          
          varcov = nu/(nu-2)*Sigma
          return(list(mean = mu,EYY = varcov + mu%*%t(mu),varcov = varcov))
          
        }else{
          
          #LEFT CENSORING
          
          return(utmvtmuvar(b = upper,mu = mu,S = Sigma,nu=nu))
        }
        
      }else{
        
        if(sum(bool2) == p){ #All upper limits are infinite?
          
          #RIGHT TRUNCATION
          
          return(ltmvtmuvar(a = lower,mu = mu,S = Sigma,nu=nu))
          
        }else{
          
          #All kind of censoring
          
          return(dtmvtmuvar(a = lower,b = upper,mu = mu,S = Sigma,nu=nu))
          
        }
        
      }
      
    }
    
  }
  
}