genPxy = function(a,b,mu,Sigma,nu){
  
  p = length(mu)
  
  if(p==2){
    
    Haa = Hab = Hba = Hbb = matrix(1,p,p)
    
    return(list(Paa = Haa, Pab = Hab, Pba = Hba, Pbb = Hbb))
    
  }
  
  Haa = Hab = Hba = Hbb = matrix(0,p,p)
  
  cdf.aa = cdf.ab = cdf.ba = cdf.bb = 1
  
  ss = sqrt(diag(as.matrix(Sigma)))
  
  aa = (a - mu)/ss
  
  bb = (b - mu)/ss
  
  R = Sigma/(ss%*%t(ss))
  
  
  for(r in 1:(p-1))
  {
    
    for(s in (r+1):p)
      
    {
      rs = c(r,s)
      
      if(p > 2)
        
      {
        
        tmp = R[-rs,rs]%*%solve(R[rs,rs])
        
        R21 = R[-rs,-rs] - R[-rs,rs]%*%solve(R[rs,rs]) %*% R[rs,-rs]
        
        
        if(all(is.finite(c(a[r],a[s]))))
          
        {
          
          mu.aa = c(tmp%*%c(aa[r],aa[s]))
          
          daa = (nu-2)/(nu+(aa[r]^2-2*R[r,s]*aa[r]*aa[s]+aa[s]^2)/(1-R[r,s]^2))
          
          cdf.aa = ifelse(p==3,pt((bb[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2)-pt((aa[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2)
                          
                          ,pmvnormt(lower = aa[-rs]-mu.aa, upper = bb[-rs]-mu.aa, sigma = R21/daa, nu=nu-2))
        }
        
        if(all(is.finite(c(a[r],b[s]))))
          
        {
          
          mu.ab = c(tmp%*%c(aa[r],bb[s]))
          
          dab = (nu-2)/(nu+(aa[r]^2-2*R[r,s]*aa[r]*bb[s]+bb[s]^2)/(1-R[r,s]^2))
          
          cdf.ab = ifelse(p==3,pt((bb[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2)-pt((aa[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2)
                          
                          ,pmvnormt(lower = aa[-rs]-mu.ab, upper = bb[-rs]-mu.ab, sigma = R21/dab, nu=nu-2))
        }
        
        if(all(is.finite(c(b[r],a[s]))))
          
        {
          
          
          mu.ba = c(tmp%*%c(bb[r],aa[s]))
          
          dba = (nu-2)/(nu+(bb[r]^2-2*R[r,s]*bb[r]*aa[s]+aa[s]^2)/(1-R[r,s]^2))
          
          cdf.ba = ifelse(p==3,pt((bb[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2)-pt((aa[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2)
                          
                          ,pmvnormt(lower = aa[-rs]-mu.ba, upper = bb[-rs]-mu.ba, sigma = R21/dba, nu=nu-2))
        }
        
        if(all(is.finite(c(b[r],b[s]))))
          
        {
          
          
          mu.bb = c(tmp%*%c(bb[r],bb[s]))
          
          dbb = (nu-2)/(nu+(bb[r]^2-2*R[r,s]*bb[r]*bb[s]+bb[s]^2)/(1-R[r,s]^2))
          
          cdf.bb = ifelse(p==3,pt((bb[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2)-pt((aa[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2)
                          
                          ,pmvnormt(lower = aa[-rs]-mu.bb, upper = bb[-rs]-mu.bb, sigma = R21/dbb, nu=nu-2))
          
        }
        
      }
      
      Haa[s,r] = Haa[r,s] = cdf.aa
      Hab[s,r] = Hba[r,s] = cdf.ab
      Hba[s,r] = Hab[r,s] = cdf.ba
      Hbb[s,r] = Hbb[r,s] = cdf.bb
      
    }
    
  }
  
  return(list(Paa = Haa, Pab = Hab, Pba = Hba, Pbb = Hbb))
}


# Upper -------------------------------------------------------------------

genPxy_upper = function(b,mu,Sigma,nu){
  
  p = length(mu)
  
  if(p==2){

    return(list(Pbb = matrix(1,p,p)))
    
  }
  
  Hbb = matrix(0,p,p)
  
  cdf.bb = 1
  
  ss = sqrt(diag(as.matrix(Sigma)))
  
  bb = (b - mu)/ss
  
  R = Sigma/(ss%*%t(ss))

  
  for(r in 1:(p-1))
    
  {
    
    for(s in (r+1):p)
      
    {
      rs = c(r,s)
      
      if(p > 2)
        
      {
        
        tmp = R[-rs,rs]%*%solve(R[rs,rs])
        
        R21 = R[-rs,-rs] - R[-rs,rs]%*%solve(R[rs,rs]) %*% R[rs,-rs]
        
        
        if(all(is.finite(b[rs])))
          
        {
          
          mu.bb = c(tmp%*%c(bb[r],bb[s]))
          
          dbb = (nu-2)/(nu+(bb[r]^2-2*R[r,s]*bb[r]*bb[s]+bb[s]^2)/(1-R[r,s]^2))
          
          cdf.bb = ifelse(p==3,pt((bb[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2),
                          
                          pmvnormt(upper = bb[-rs]-mu.bb, sigma = R21/dbb, nu=nu-2))
          
        }
        
      }
      
      Hbb[s,r] = Hbb[r,s] = cdf.bb
      
    }
    
  }
  
  return(list(Pbb = Hbb))
}


# Finite ------------------------------------------------------------------

genPxy_finite = function(a,b,mu,Sigma,nu){
  
  p = length(mu)
  
  if(p==2){
    
    Haa = Hab = Hba = Hbb = matrix(1,p,p)
    
    return(list(Paa = Haa, Pab = Hab, Pba = Hba, Pbb = Hbb))
    
  }
  
  Haa = Hab = Hba = Hbb = matrix(0,p,p)
  
  cdf.aa = cdf.ab = cdf.ba = cdf.bb = 1
  
  ss = sqrt(diag(as.matrix(Sigma)))
  
  aa = (a - mu)/ss
  
  bb = (b - mu)/ss
  
  R = Sigma/(ss%*%t(ss))
  
  
  for(r in 1:(p-1))
    
  {
    
    for(s in (r+1):p)
      
    {
      rs = c(r,s)
      
      if(p > 2)
        
      {
        
        tmp = R[-rs,rs]%*%solve(R[rs,rs])
        
        R21 = R[-rs,-rs] - R[-rs,rs]%*%solve(R[rs,rs]) %*% R[rs,-rs]

        
        mu.aa = c(tmp%*%c(aa[r],aa[s]))
        
        daa = (nu-2)/(nu+(aa[r]^2-2*R[r,s]*aa[r]*aa[s]+aa[s]^2)/(1-R[r,s]^2))
        
        cdf.aa = ifelse(p==3,pt((bb[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2)-pt((aa[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2)
                        
                        ,pmvnormt(lower = aa[-rs]-mu.aa, upper = bb[-rs]-mu.aa, sigma = R21/daa, nu=nu-2))
        
        
        mu.ab = c(tmp%*%c(aa[r],bb[s]))
        
        dab = (nu-2)/(nu+(aa[r]^2-2*R[r,s]*aa[r]*bb[s]+bb[s]^2)/(1-R[r,s]^2))
        
        cdf.ab = ifelse(p==3,pt((bb[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2)-pt((aa[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2)
                        
                        ,pmvnormt(lower = aa[-rs]-mu.ab, upper = bb[-rs]-mu.ab, sigma = R21/dab, nu=nu-2))
        
        
        
        mu.ba = c(tmp%*%c(bb[r],aa[s]))
        
        dba = (nu-2)/(nu+(bb[r]^2-2*R[r,s]*bb[r]*aa[s]+aa[s]^2)/(1-R[r,s]^2))
        
        cdf.ba = ifelse(p==3,pt((bb[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2)-pt((aa[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2)
                        
                        ,pmvnormt(lower = aa[-rs]-mu.ba, upper = bb[-rs]-mu.ba, sigma = R21/dba, nu=nu-2))
        
        
        
        mu.bb = c(tmp%*%c(bb[r],bb[s]))
        
        dbb = (nu-2)/(nu+(bb[r]^2-2*R[r,s]*bb[r]*bb[s]+bb[s]^2)/(1-R[r,s]^2))
        
        cdf.bb = ifelse(p==3,pt((bb[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2)-pt((aa[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2)
                        
                        ,pmvnormt(lower = aa[-rs]-mu.bb, upper = bb[-rs]-mu.bb, sigma = R21/dbb, nu=nu-2))
        
      }
      
      Haa[s,r] = Haa[r,s] = cdf.aa
      Hab[s,r] = Hba[r,s] = cdf.ab
      Hba[s,r] = Hab[r,s] = cdf.ba
      Hbb[s,r] = Hbb[r,s] = cdf.bb
      
    }
    
  }
  
  return(list(Paa = Haa, Pab = Hab, Pba = Hba, Pbb = Hbb))
}

