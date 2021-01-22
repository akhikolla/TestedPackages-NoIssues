###################################################################################
###################################################################################
#AUXILIAR FUNCTIONS
###################################################################################
###################################################################################

qfun = function(a,b,Sigma)
{
  n = length(a)
  s = sqrt(diag(as.matrix(Sigma)))
  if(n==1){
    qa = dnorm(a/s)/s
    qb = dnorm(b/s)/s
    return(list(qa = qa,qb = qb))
  }
  if(n==2){
    qa = qb = rep(0,n)
    for(i in 1:n){
      if(a[i] != -Inf){
        qa[i] = dnorm(x = a[i],mean = 0,sd = s[i])*pnorm2(lower = a[-i],upper = b[-i],mean = c(Sigma[-i,i]/Sigma[i,i]*a[i]),sd = sqrt(c(Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i]))))
      }
      if(b[i] != Inf){
        qb[i] = dnorm(x = b[i],mean = 0,sd = s[i])*pnorm2(lower = a[-i],upper = b[-i],mean = c(Sigma[-i,i]/Sigma[i,i]*b[i]),sd = sqrt(c(Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i]))))
      }
    }
    return(list(qa = qa,qb = qb))
  }
  qa = qb = rep(0,n)
  for(i in 1:n){
    if(a[i] != -Inf){
      qa[i] = log2prod0(dnorm(x = a[i],mean = 0,sd = s[i]),
                        pmvnormt(lower = a[-i],upper = b[-i],mean = Sigma[-i,i]/Sigma[i,i]*a[i],sigma = Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i]),uselog2 = TRUE))
    }
    if(b[i] != Inf){
      qb[i] = log2prod0(dnorm(x = b[i],mean = 0,sd = s[i]),
                        pmvnormt(lower = a[-i],upper = b[-i],mean = Sigma[-i,i]/Sigma[i,i]*b[i],sigma = Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i]),uselog2 = TRUE))
    }
  }
  return(list(qa = qa,qb = qb))
}

###################################################################################
###################################################################################

qfun_b = function(b1,Sigma)
{
  p = length(b1)
  s = sqrt(diag(as.matrix(Sigma)))
  if(p==1){
    qb = dnorm(b1/s)/s
    return(qb)
  }
  if(p==2){
    qb = rep(0,p)
    for(i in 1:p){
      if(b1[i] != Inf){
        qb[i] = dnorm(x = b1[i],mean = 0,sd = s[i])*pnorm2(upper = b1[-i],mean = c(Sigma[-i,i]/Sigma[i,i]*b1[i]),sd = sqrt(c(Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i]))))
      }
    }
    return(qb)
  }
  qb = rep(0,p)
  for(i in 1:p){
    if(b1[i] != Inf){
      qb[i] = log2prod0(dnorm(x = b1[i],mean = 0,sd = s[i]),
                        pmvnormt(upper = b1[-i],mean = Sigma[-i,i]/Sigma[i,i]*b1[i],sigma = Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i]),uselog2 = TRUE))
    }
  }
  return(qb)
}

###################################################################################
###################################################################################

qfun.noinf = function(a,b,Sigma)
{
  n = length(a)
  s = sqrt(diag(as.matrix(Sigma)))
  if(n==1){
    qa = dnorm(a/s)/s
    qb = dnorm(b/s)/s
    return(list(qa = qa,qb = qb))
  }
  if(n==2){
    qa = qb = rep(0,n)
    for(i in 1:n){
      qa[i] = dnorm(x = a[i],mean = 0,sd = s[i])*pnorm2(lower = a[-i],upper = b[-i],mean = c(Sigma[-i,i]/Sigma[i,i]*a[i]),sd = sqrt(c(Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i]))))
      qb[i] = dnorm(x = b[i],mean = 0,sd = s[i])*pnorm2(lower = a[-i],upper = b[-i],mean = c(Sigma[-i,i]/Sigma[i,i]*b[i]),sd = sqrt(c(Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i]))))
    }
    return(list(qa = qa,qb = qb))
  }
  qa = qb = rep(0,n)
  for(i in 1:n){
    qa[i] = log2prod0(dnorm(x = a[i],mean = 0,sd = s[i]),
                      pmvnormt(lower = a[-i],upper = b[-i],mean = Sigma[-i,i]/Sigma[i,i]*a[i],sigma = Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i]),uselog2 = TRUE))
    qb[i] = log2prod0(dnorm(x = b[i],mean = 0,sd = s[i]),
                      pmvnormt(lower = a[-i],upper = b[-i],mean = Sigma[-i,i]/Sigma[i,i]*b[i],sigma = Sigma[-i,-i] - (Sigma[-i,i]/Sigma[i,i])%*%t(Sigma[i,-i]),uselog2 = TRUE))
  }
  return(list(qa = qa,qb = qb))
}
