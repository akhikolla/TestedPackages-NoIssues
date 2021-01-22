###############################################################################################
###############################################################################################
#This function is optimized for the case when all intervalar lower/upper censoring limits
#are finite
###############################################################################################
###############################################################################################

Kan.IC = function(a,b,mu,Sigma){
  n = length(mu)
  s = sqrt(diag(Sigma))
  seqq = seq_len(n)
  a1 = a-mu
  b1 = b-mu

  #####
  Sigma = sym.matrix(Sigma)
  #####

  logp = pmvnormt(lower = a,upper = b,mean = mu,sigma = Sigma,uselog2 = TRUE)
  
  prob = 2^logp
  # if(prob > 1e-50){
  #   #no problems, so we run the Rcpp model
  #   #print("no problems")
  #   print("Rcpp")
  #   return(RcppmeanvarN_ab(a,b,mu,Sigma,prob))
  # }
  
  if(prob < 1e-250){
    #print("corrector")
    #print("LRIC.Kan corrector applied \n")
    return(corrector(a,b,mu,Sigma,bw=36))
  }
  
  run = Rcppqfun_ab(a1,b1,Sigma)
  qa = run$qa
  qb = run$qb
  q = qa-qb
  muY = mu+ Sigma%*%log2ratio(q,logp)

  if(
    #max(abs(muY)) > 10*max(abs(c(a,b)[is.finite(c(a,b))]))
    any(b < mu - 10*sqrt(diag(Sigma))) |
    any(a > mu + 10*sqrt(diag(Sigma))) |
    any(muY < a | muY > b)){
    #print("IC.Kan mean corrector applied 2 \n")
    #print("corrector")
    return(corrector(a,b,mu,Sigma,bw=36))
  }

  D = matrix(0,n,n)
  for(i in seqq){
    D[i,i] = a[i]*qa[i]
    D[i,i] = D[i,i]-b[i]*qb[i]
    RR = Sigma[-i,-i]-Sigma[-i,i]%*%t(Sigma[i,-i])/Sigma[i,i]
    ma = mu[-i]+Sigma[-i,i]/Sigma[i,i]*a1[i]
    run1 = Rcppqfun_ab(a[-i]-ma,b[-i]-ma,RR)
    qa1 = run1$qa
    qb1 = run1$qb
    wa = qa[i]*ma+dnorm(x = a[i],mean = mu[i],sd = s[i])*RR%*%(qa1-qb1)
    mb = mu[-i]+Sigma[-i,i]/Sigma[i,i]*b1[i]
    run2 = Rcppqfun_ab(a[-i]-mb,b[-i]-mb,RR)
    qa2 = run2$qa
    qb2 = run2$qb
    wb = qb[i]*mb + dnorm(x = b[i],mean = mu[i],sd = s[i])*RR%*%(qa2-qb2)
    D[i,-i] = wa-wb
  }
  varY = Sigma + Sigma%*%log2ratio(D - q%*%t(muY),logp)
  varY = (varY + t(varY))/2
  EYY = varY+muY%*%t(muY)

  bool = diag(varY) < 0
  if(sum(bool)>0){
    
    #print("corrector")
    
    #print("negative variance found")
    out = corrector(a,b,mu,Sigma,bw=36)
    out$mean = muY
    out$EYY = out$varcov + out$mean%*%t(out$mean)
    return(out)
  }

  return(list(mean = muY,EYY = EYY,varcov = varY))
}

###############################################################################################
###############################################################################################
#This function is optimized for the case when it DOES exist infinite values in the lower/upper
#truncation limits
###############################################################################################
###############################################################################################

Kan.LRIC = function(a,b,mu,Sigma){
  n = length(mu)
  s = sqrt(diag(Sigma))
  seqq = seq_len(n)
  a1 = a-mu
  b1 = b-mu

  #####
  Sigma = sym.matrix(Sigma)
  #####

  logp = pmvnormt(lower = a,upper = b,mean = mu,sigma = Sigma,uselog2 = TRUE)
  prob = 2^logp
  # if(prob > 1e-50){
  #   #no problems, so we run the Rcpp model
  #   print("Rcpp")
  #   return(RcppmeanvarN(a,b,mu,Sigma,prob))
  # }
  
  if(prob < 1e-250){
    #print("corrector")
    #print("LRIC.Kan corrector applied \n")
    return(corrector(a,b,mu,Sigma,bw=36))
  }
  
  run = Rcppqfun(a1,b1,Sigma)
  qa = run$qa
  qb = run$qb
  q = qa-qb
  muY = mu+ Sigma%*%log2ratio(q,logp)

  if(any(b < mu - 10*sqrt(diag(Sigma))) |
     any(a > mu + 10*sqrt(diag(Sigma))) | any(muY < a | muY > b)){
    #print("corrector")
    return(corrector(a,b,mu,Sigma,bw=36))
  }

  D = matrix(0,n,n)
  for(i in seqq){
    if(a[i] != -Inf){
      D[i,i] = a[i]*qa[i]
    }
    if(b[i] != Inf){
      D[i,i] = D[i,i]-b[i]*qb[i]
    }
    RR = Sigma[-i,-i]-Sigma[-i,i]%*%t(Sigma[i,-i])/Sigma[i,i]
    if(a[i] == -Inf){
      wa = matrix(0,n-1,1)
    }else
    {
      ma = mu[-i]+Sigma[-i,i]/Sigma[i,i]*a1[i]
      run1 = Rcppqfun(a[-i]-ma,b[-i]-ma,RR)
      qa1 = run1$qa
      qb1 = run1$qb
      wa = qa[i]*ma+dnorm(x = a[i],mean = mu[i],sd = s[i])*RR%*%(qa1-qb1)
    }
    if(b[i] == Inf){
      wb = matrix(0,n-1,1)
    }else
    {
      mb = mu[-i]+Sigma[-i,i]/Sigma[i,i]*b1[i]
      run2 = Rcppqfun(a[-i]-mb,b[-i]-mb,RR)
      qa2 = run2$qa
      qb2 = run2$qb
      wb = qb[i]*mb + dnorm(x = b[i],mean = mu[i],sd = s[i])*RR%*%(qa2-qb2)
    }
    D[i,-i] = wa-wb
  }
  varY = Sigma + Sigma%*%log2ratio(D - q%*%t(muY),logp)
  varY = (varY + t(varY))/2
  EYY = varY+muY%*%t(muY)

  #Validating positive variances
  bool = diag(varY) < 0
  if(sum(bool)>0){
    #print("corrector")
    #print("negative variance found")
    out = corrector(a,b,mu,Sigma,bw=36)
    out$mean = muY
    out$EYY = out$varcov + out$mean%*%t(out$mean)
    return(out)
  }

  return(list(mean = muY,EYY = EYY,varcov = varY))
}

###############################################################################################
###############################################################################################
#right censoring
###############################################################################################
###############################################################################################

Kan.RC = function(b,mu,Sigma){
  n = length(mu)
  s = sqrt(diag(Sigma))
  seqq = seq_len(n)
  b1 = b-mu

  #####
  Sigma = sym.matrix(Sigma)
  #####

  logp = pmvnormt(upper = b,mean = mu,sigma = Sigma,uselog2 = TRUE)
  
  prob = 2^logp
  # if(prob > 1e-50){
  #   #no problems, so we run the Rcpp model
  #   print("Rcpp")
  #   return(RcppmeanvarN_b(b,mu,Sigma,prob))
  # }
  
  if(prob < 1e-250){
    #print("corrector")
    #print("LRIC.Kan corrector applied \n")
    return(corrector(upper = b,mu = mu,Sigma = Sigma,bw=36))
  }
  
  qb = Rcppqfun_b(b1,Sigma)
  muY = mu - Sigma%*%log2ratio(qb,logp)

#max(abs(muY))> 10*max(abs(b[is.finite(b)]))
  
  if(any(b < mu - 10*sqrt(diag(Sigma))) | any(muY > b)){
    #print("corrector")
    return(corrector(upper = b,mu = mu,Sigma = Sigma,bw=36))
  }

  D = matrix(0,n,n)
  for(i in seqq){
    D[i,i] = D[i,i]-b[i]*qb[i]
    RR = Sigma[-i,-i]-Sigma[-i,i]%*%t(Sigma[i,-i])/Sigma[i,i]
    mb = mu[-i]+Sigma[-i,i]/Sigma[i,i]*b1[i]
    qb2 = Rcppqfun_b(b[-i]-mb,RR)
    wb = qb[i]*mb - dnorm(x = b[i],mean = mu[i],sd = s[i])*RR%*%qb2
    D[i,-i] = -wb
  }
  varY = Sigma + Sigma%*%log2ratio(D + qb%*%t(muY),logp)
  varY = (varY + t(varY))/2
  EYY = varY+muY%*%t(muY)

  #Validating positive variances
  bool = diag(varY) < 0
  if(sum(bool)>0){
    #print("corrector")
    #print("negative variance found")
    out = corrector(upper = b,mu = mu,Sigma = Sigma,bw=36)
    out$mean = muY
    out$EYY = out$varcov + out$mean%*%t(out$mean)
    return(out)
  }

  return(list(mean = muY,EYY = EYY,varcov = varY))
}

# ########################
# #TESTING
# ########################
#
# p = 4
# mu  = c(matrix(c(1:p/10)))
# s  = 2*matrix(rnorm(p^2),p,p)
# Sigma = S = round(s%*%t(s)/10,1)
# lambda = seq(-1,2,length.out = p)
# tau = 1
#
# a = rep(-Inf,p)
# b = mu+2
#
# Kan.R(b,mu,Sigma)
# Vaida(b,mu,Sigma)
#
# compare <- microbenchmark(Kan.R(b,mu,Sigma),
#                           Vaida(b,mu,Sigma),
#                           times = 100)
# autoplot(compare)
#
# ########################
#
# a = mu-1
# b = mu+2
# a[2] = -Inf
# b[3] = Inf
#
# Kan.LRIC(a,b,mu,Sigma)
# meanvarN(a,b,mu,Sigma)
