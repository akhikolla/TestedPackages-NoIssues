# R commands: calculation of the first two moments of the TMVT distribution

meanvarT.Lin.LRIC = function(a,b,mu,S,nu,omega = FALSE){
  ss = sqrt(diag(S))
  as = (a-mu)/ss
  bs = (b-mu)/ss
  R = S/(ss%*%t(ss))
  M=TT.moment.LRIC(R = R,nu, lower=as, upper=bs)
  M$varcov = M$EYY -  M$mean%*%t(M$mean)
  M$mean = ss*M$mean + mu
  M$varcov = diag(ss)%*%M$varcov%*%diag(ss)
  M$EYY = M$varcov + M$mean%*%t(M$mean)
  if(omega){
    return(M)
  }else{
    return(list(mean = M$mean, EYY = M$EYY, varcov = M$varcov))
  }
}


meanvarT.Lin.IC = function(a,b,mu,S,nu,omega = FALSE){
  ss = sqrt(diag(S))
  as = (a-mu)/ss
  bs = (b-mu)/ss
  R = S/(ss%*%t(ss))
  M=TT.moment.IC(R = R,nu, lower=as, upper=bs)
  M$varcov = M$EYY -  M$mean%*%t(M$mean)
  M$mean = ss*M$mean + mu
  M$varcov = diag(ss)%*%M$varcov%*%diag(ss)
  M$EYY = M$varcov + M$mean%*%t(M$mean)
  if(omega){
    return(M)
  }else{
    return(list(mean = M$mean, EYY = M$EYY, varcov = M$varcov))
  }
}


meanvarT.Lin.RC = function(b,mu,S,nu, omega = FALSE){
  ss = sqrt(diag(S))
  bs = (b-mu)/ss
  R = S/(ss%*%t(ss))
  M=TT.moment.RC(R = R,nu,upper=bs)
  
  M = TTmoment::TT.moment(R,nu,lower = c(-Inf,-Inf),upper = bs)
  M$mean = M$EX
  M$EYY = M$EXX
  
  M = TT.moment.LRIC(R = R,nu = nu,lower = c(-Inf,-Inf),upper = bs)
  
  
  M$varcov = M$EYY -  M$mean%*%t(M$mean)
  M$mean = ss*M$mean + mu
  M$varcov = diag(ss)%*%M$varcov%*%diag(ss)
  M$EYY = M$varcov + M$mean%*%t(M$mean)
  if(omega){
    return(M)
  }else{
    return(list(mean = M$mean, EYY = M$EYY, varcov = M$varcov))
  }
}

meanvarT.Lin.LC = function(a,mu,S,nu, omega = FALSE){
  out = meanvarT.Lin.RC(-a,-mu,S,nu,omega)
  out$mean = -out$mean
  return(out)
}


# All ---------------------------------------------------------------------



TT.moment.LRIC = function(R=diag(length(lower)), nu=5, lower=rep(-Inf, nrow(R)), upper=rep(Inf, nrow(R)))
  
{
  
  a = lower; b = upper
  p = length(a)
  
  al0 = pmvnormt(lower = a, upper = b, sigma = R, nu = nu,uselog2 = TRUE)
  
  ### pdf & cdf
  
  la1 = (nu-2)/nu; la2 = (nu-4)/nu
  
  da = (nu-1)/(nu+a^2); db = (nu-1)/(nu+b^2)
  
  f1a = sqrt(la1)*dt(sqrt(la1)*a,df=nu-2)
  
  f1b = sqrt(la1)*dt(sqrt(la1)*b,df=nu-2)
  
  f2 = matrix(NA, p, p)
  
  G1a = G1b = rep(0, p)
  
  flag.a = is.finite(a)
  flag.b = is.finite(b)
  
  for(r in 1:p)
    
  {
    
    temp = R[-r,r]
    
    S1 = R[-r,-r] - temp %*% t(R[r,-r])
    
    if(flag.a[r]){
      
      mua = temp * a[r]; low = a[-r]-mua; upp = b[-r]-mua
      
      G1a[r] = ifelse(p==2,pt(upp/sqrt(S1/da[r]),df=nu-1)-pt(low/sqrt(S1/da[r]),df=nu-1)
                      
                      ,pmvnormt(lower = low, upper = upp, sigma = S1/da[r], nu = nu-1))
    }
    
    if(flag.b[r]){
      
      mub = temp * b[r]; low = a[-r]-mub; upp = b[-r]-mub
      
      G1b[r] = ifelse(p==2,pt(upp/sqrt(S1/db[r]),df=nu-1)-pt(low/sqrt(S1/db[r]),df=nu-1)
                      
                      ,pmvnormt(lower = low, upper = upp, sigma = S1/db[r], nu = nu-1))
    }
    
  }
  
  qa = f1a*G1a; qb = f1b*G1b
  
  EX = R %*% log2ratio(qa-qb,al0) / la1
  
  if(nu>4){
    
    H = matrix(0,p,p)
    
    cdf.aa = cdf.ab = cdf.ba = cdf.bb = 0
    
    for(r in 1:(p-1))
      
    {
      
      for(s in (r+1):p)
        
      {
        rs = c(r,s)
        
        # pdf.aa = bivT(c(a[r],a[s]),S=R[rs,rs]/la2,nu=nu-4)
        # 
        # pdf.ab = bivT(c(a[r],b[s]),S=R[rs,rs]/la2,nu=nu-4)
        # 
        # pdf.ba = bivT(c(b[r],a[s]),S=R[rs,rs]/la2,nu=nu-4)
        # 
        # pdf.bb = bivT(c(b[r],b[s]),S=R[rs,rs]/la2,nu=nu-4)
        
        
        pdf.aa = dt2d(c(a[r],a[s])*sqrt(la2),rho = R[r,s],nu = nu-4)*la2
        
        pdf.ab = dt2d(c(a[r],b[s])*sqrt(la2),rho = R[r,s],nu = nu-4)*la2
        
        pdf.ba = dt2d(c(b[r],a[s])*sqrt(la2),rho = R[r,s],nu = nu-4)*la2
        
        pdf.bb = dt2d(c(b[r],b[s])*sqrt(la2),rho = R[r,s],nu = nu-4)*la2
        
        
        if(p==2){cdf.aa = cdf.ab = cdf.ba = cdf.bb = 1}
        
        if(p > 2)
          
        {
          
          tmp = R[-rs,rs]%*%solve(R[rs,rs])
          
          R21 = R[-rs,-rs] - R[-rs,rs]%*%solve(R[rs,rs]) %*% R[rs,-rs]
          
          
          if(pdf.aa != 0)
            
          {
  
            mu.aa = c(tmp%*%c(a[r],a[s]))
            
            daa = (nu-2)/(nu+(a[r]^2-2*R[r,s]*a[r]*a[s]+a[s]^2)/(1-R[r,s]^2))
            
            cdf.aa = ifelse(p==3,pt((b[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2)-pt((a[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2)
                            
                            ,pmvnormt(lower = a[-rs]-mu.aa, upper = b[-rs]-mu.aa, sigma = R21/daa, nu=nu-2))
          }
          
          if(pdf.ab != 0)
            
          {
            
            mu.ab = c(tmp%*%c(a[r],b[s]))
            
            dab = (nu-2)/(nu+(a[r]^2-2*R[r,s]*a[r]*b[s]+b[s]^2)/(1-R[r,s]^2))
            
            cdf.ab = ifelse(p==3,pt((b[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2)-pt((a[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2)
                            
                            ,pmvnormt(lower = a[-rs]-mu.ab, upper = b[-rs]-mu.ab, sigma = R21/dab, nu=nu-2))
          }
          
          if(pdf.ba != 0)
            
          {
            
            
            mu.ba = c(tmp%*%c(b[r],a[s]))
            
            dba = (nu-2)/(nu+(b[r]^2-2*R[r,s]*b[r]*a[s]+a[s]^2)/(1-R[r,s]^2))
            
            cdf.ba = ifelse(p==3,pt((b[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2)-pt((a[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2)
                            
                            ,pmvnormt(lower = a[-rs]-mu.ba, upper = b[-rs]-mu.ba, sigma = R21/dba, nu=nu-2))
          }
          
          if(pdf.bb != 0)
            
          {
            
            
            mu.bb = c(tmp%*%c(b[r],b[s]))
            
            dbb = (nu-2)/(nu+(b[r]^2-2*R[r,s]*b[r]*b[s]+b[s]^2)/(1-R[r,s]^2))
            
            cdf.bb = ifelse(p==3,pt((b[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2)-pt((a[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2)
                            
                            ,pmvnormt(lower = a[-rs]-mu.bb, upper = b[-rs]-mu.bb, sigma = R21/dbb, nu=nu-2))
            
          }
          
        }
        
        H[r,s] = H[s,r] = pdf.aa*cdf.aa - pdf.ab*cdf.ab - pdf.ba*cdf.ba + pdf.bb*cdf.bb
        
      }
      
    }
    
    H = H / la2
    
    D = matrix(0,p,p)
    
    al1 = pmvnormt(lower = a, upper = b, sigma = R/la1, nu=nu-2,uselog2 = TRUE)
    
    a[is.infinite(a)] = 0
    b[is.infinite(b)] = 0
    
    diag(D) = a * qa - b * qb - diag(R%*%H)
    
    #EXX = (2^al1 * R + R %*% (H + D) %*% R) / 2^al0 / la1
    #EXX = log2ratio(2^al1 * R + R %*% (H + D) %*% R, al0) / la1
    #EXX = log2ratio(2^al1 * R, al0) / la1 + log2ratio(R %*% (H + D) %*% R, al0) / la1
    
    ratio0 = 2^(al1-al0)
    
    #EXX = (al1 * R + R %*% (H + D) %*% R) / al0 / la1
    EXX = (ratio0 * R + log2ratio(R %*% (H + D) %*% R, al0)) / la1
    omega0 = nu/(nu-2)*ratio0
    
  } else {
    EXX=matrix(NA,p,p)
    omega0 = NA
    cat('Warning message:','\n')
    cat('It only works for degrees of freedom larger than 4!')
  }
  
  return(list(mean=EX,EYY=EXX,omega = omega0))
  
}

# Intervalar --------------------------------------------------------------



TT.moment.IC = function(R=diag(length(lower)), nu=5, lower=rep(-Inf, nrow(R)), upper=rep(Inf, nrow(R)))
  
{
  
  a = lower; b = upper
  p = length(a)
  
  al0 = pmvnormt(lower = a, upper = b, sigma = R, nu = nu,uselog2 = TRUE)
  
  ### pdf & cdf
  
  la1 = (nu-2)/nu; la2 = (nu-4)/nu
  
  da = (nu-1)/(nu+a^2); db = (nu-1)/(nu+b^2)
  
  f1a = sqrt(la1)*dt(sqrt(la1)*a,df=nu-2)
  
  f1b = sqrt(la1)*dt(sqrt(la1)*b,df=nu-2)
  
  f2 = matrix(NA, p, p)
  
  G1a = G1b = rep(0, p)
  
  
  for(r in 1:p)
    
  {
    
    temp = R[-r,r]
    
    S1 = R[-r,-r] - temp %*% t(R[r,-r])
    
    mua = temp * a[r]; low = a[-r]-mua; upp = b[-r]-mua
    
    G1a[r] = ifelse(p==2,pt(upp/sqrt(S1/da[r]),df=nu-1)-pt(low/sqrt(S1/da[r]),df=nu-1)
                    
                    ,pmvnormt(lower = low, upper = upp, sigma = S1/da[r], nu = nu-1))
    
    mub = temp * b[r]; low = a[-r]-mub; upp = b[-r]-mub
    
    G1b[r] = ifelse(p==2,pt(upp/sqrt(S1/db[r]),df=nu-1)-pt(low/sqrt(S1/db[r]),df=nu-1)
                    
                    ,pmvnormt(lower = low, upper = upp, sigma = S1/db[r], nu = nu-1))
    
  }
  
  qa = f1a*G1a; qb = f1b*G1b
  
  EX = R %*% log2ratio(qa-qb,al0) / la1
  
  if(nu>4){
    
    H = matrix(0,p,p)
    
    for(r in 1:(p-1))
      
    {
      
      for(s in (r+1):p)
        
      {
        
        rs = c(r,s)
        
        
        # pdf.aa = bivT(c(a[r],a[s]),S=R[rs,rs]/la2,nu=nu-4)
        # 
        # pdf.ab = bivT(c(a[r],b[s]),S=R[rs,rs]/la2,nu=nu-4)
        # 
        # pdf.ba = bivT(c(b[r],a[s]),S=R[rs,rs]/la2,nu=nu-4)
        # 
        # pdf.bb = bivT(c(b[r],b[s]),S=R[rs,rs]/la2,nu=nu-4)
        
        
        pdf.aa = dt2d(c(a[r],a[s])*sqrt(la2),rho = R[r,s],nu = nu-4)*la2
        
        pdf.ab = dt2d(c(a[r],b[s])*sqrt(la2),rho = R[r,s],nu = nu-4)*la2
        
        pdf.ba = dt2d(c(b[r],a[s])*sqrt(la2),rho = R[r,s],nu = nu-4)*la2
        
        pdf.bb = dt2d(c(b[r],b[s])*sqrt(la2),rho = R[r,s],nu = nu-4)*la2
        
        if(p==2){cdf.aa=cdf.ab=cdf.ba=cdf.bb=1}
        
        if(p>2)
          
        {
          
          tmp = R[-rs,rs]%*%solve(R[rs,rs])
          
          mu.aa = c(tmp%*%c(a[r],a[s]))
          
          mu.ab = c(tmp%*%c(a[r],b[s]))
          
          mu.ba = c(tmp%*%c(b[r],a[s]))
          
          mu.bb = c(tmp%*%c(b[r],b[s]))
          
          daa = (nu-2)/(nu+(a[r]^2-2*R[r,s]*a[r]*a[s]+a[s]^2)/(1-R[r,s]^2))
          
          dab = (nu-2)/(nu+(a[r]^2-2*R[r,s]*a[r]*b[s]+b[s]^2)/(1-R[r,s]^2))
          
          dba = (nu-2)/(nu+(b[r]^2-2*R[r,s]*b[r]*a[s]+a[s]^2)/(1-R[r,s]^2))
          
          dbb = (nu-2)/(nu+(b[r]^2-2*R[r,s]*b[r]*b[s]+b[s]^2)/(1-R[r,s]^2))
          
          R21 = R[-rs,-rs] - R[-rs,rs]%*%solve(R[rs,rs]) %*% R[rs,-rs]
          
          cdf.aa = ifelse(p==3,pt((b[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2)-pt((a[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2)
                          
                          ,pmvnormt(lower = a[-rs]-mu.aa, upper = b[-rs]-mu.aa, sigma = R21/daa, nu=nu-2))
          
          cdf.ab = ifelse(p==3,pt((b[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2)-pt((a[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2)
                          
                          ,pmvnormt(lower = a[-rs]-mu.ab, upper = b[-rs]-mu.ab, sigma = R21/dab, nu=nu-2))
          
          cdf.ba = ifelse(p==3,pt((b[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2)-pt((a[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2)
                          
                          ,pmvnormt(lower = a[-rs]-mu.ba, upper = b[-rs]-mu.ba, sigma = R21/dba, nu=nu-2))
          
          cdf.bb = ifelse(p==3,pt((b[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2)-pt((a[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2)
                          
                          ,pmvnormt(lower = a[-rs]-mu.bb, upper = b[-rs]-mu.bb, sigma = R21/dbb, nu=nu-2))
          
        }
        
        H[r,s] = H[s,r] = pdf.aa*cdf.aa - pdf.ab*cdf.ab - pdf.ba*cdf.ba + pdf.bb*cdf.bb
        
      }
      
    }
    
    H = H / la2
    
    D = matrix(0,p,p)
    
    diag(D) = a * qa - b * qb - diag(R%*%H)
    
    al1 = pmvnormt(lower = a, upper = b, sigma = R/la1, nu=nu-2,uselog2 = TRUE)
    
    #EXX = (2^al1 * R + R %*% (H + D) %*% R) / 2^al0 / la1
    #EXX = log2ratio(2^al1 * R + R %*% (H + D) %*% R, al0) / la1
    #EXX = log2ratio(2^al1 * R, al0) / la1 + log2ratio(R %*% (H + D) %*% R, al0) / la1
    ratio0 = 2^(al1-al0)
    EXX = (ratio0 * R + log2ratio(R %*% (H + D) %*% R, al0)) / la1
    
    omega0 = nu/(nu-2)*ratio0
    
  } else {
    EXX=matrix(NA,p,p)
    omega0 = NA
    cat('Warning message:','\n')
    cat('It only works for degrees of freedom larger than 4!')
  }
  
  return(list(mean=EX,EYY=EXX,omega = omega0))
  
}


# Right censoring ---------------------------------------------------------


TT.moment.RC = function(R=diag(length(upper)), nu=5, upper=rep(Inf, nrow(R)))
  
{
  
  b = upper
  p = length(b)
  
  al0 = pmvnormt(upper = b, sigma = R, nu = nu,uselog2 = TRUE)
  
  ### pdf & cdf
  
  la1 = (nu-2)/nu; la2 = (nu-4)/nu
  
  db = (nu-1)/(nu+b^2)
  
  f1b = sqrt(la1)*dt(sqrt(la1)*b,df=nu-2)
  
  G1b = rep(0, p)
  
  flag.b = is.finite(b)
  
  for(r in 1:p)
    
  {
    
    temp = R[-r,r]
    
    S1 = R[-r,-r] - temp %*% t(R[r,-r])
    
    
    if(flag.b[r]){
      
      mub = temp * b[r]; upp = b[-r]-mub
      
      G1b[r] = ifelse(p==2,pt(upp/sqrt(S1/db[r]),df=nu-1),
                      pmvnormt(upper = upp, sigma = S1/db[r], nu = nu-1))
    }
    
  }
  
  qb = f1b*G1b
  
  EX = - log2ratio(R%*%qb,al0) / la1
  
  if(nu>4){
    
    H = matrix(0,p,p)
    
    cdf.bb = 0
    
    for(r in 1:(p-1))
      
    {
      
      for(s in (r+1):p)
        
      {
        rs = c(r,s)
        
        # pdf.aa = bivT(c(a[r],a[s]),S=R[rs,rs]/la2,nu=nu-4)
        # 
        # pdf.ab = bivT(c(a[r],b[s]),S=R[rs,rs]/la2,nu=nu-4)
        # 
        # pdf.ba = bivT(c(b[r],a[s]),S=R[rs,rs]/la2,nu=nu-4)
        # 
        # pdf.bb = bivT(c(b[r],b[s]),S=R[rs,rs]/la2,nu=nu-4)
        
        
        pdf.bb = dt2d(c(b[r],b[s])*
                        sqrt(la2),rho = R[r,s],nu = nu-4)*la2
        
        
        if(p==2){cdf.bb = 1}
        
        if(p > 2)
          
        {
          
          tmp = R[-rs,rs]%*%solve(R[rs,rs])
          
          R21 = R[-rs,-rs] - R[-rs,rs]%*%solve(R[rs,rs]) %*% R[rs,-rs]
          
          
          if(pdf.bb != 0)
            
          {
            
            
            mu.bb = c(tmp%*%c(b[r],b[s]))
            
            dbb = (nu-2)/(nu+(b[r]^2-2*R[r,s]*b[r]*b[s]+b[s]^2)/(1-R[r,s]^2))
            
            cdf.bb = ifelse(p==3,pt((b[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2),
                            pmvnormt(upper = b[-rs]-mu.bb, sigma = R21/dbb, nu=nu-2))
            
          }
          
        }
        
        H[r,s] = H[s,r] = pdf.bb*cdf.bb
        
      }
      
    }
    
    H = H / la2
    
    D = matrix(0,p,p)
    
    al1 = pmvnormt(upper = b, sigma = R/la1, nu=nu-2,uselog2 = TRUE)
    
    b[is.infinite(b)] = 0
    
    diag(D) = - b * qb - diag(R%*%H)
    
    
    ratio0 = 2^(al1-al0)
    EXX = (ratio0 * R + log2ratio(R %*% (H + D) %*% R, al0)) / la1
    omega0 = nu/(nu-2)*ratio0
    
    EXX - EX%*%t(EX)
    
  } else {
    EXX=matrix(NA,p,p)
    omega0 = NA
    cat('Warning message:','\n')
    cat('It only works for degrees of freedom larger than 4!')
  }
  
  return(list(mean=EX,EYY=EXX,omega = omega0))
  
}
