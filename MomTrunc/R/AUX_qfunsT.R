# %
# %   This function computes two vectors of integrals, qa and qb.
# %   When qa(i) and qb(i) do not exist but qa(i)-qb(i) exists,
# %   the answer is stored in qa(i).  qa(i) and qb(i) are
# %   defined as
# %   qa(i) = c_{alpha_i}/sigma_i*L_{n-1,nu-1}(a(i),b(i),\tilde\mu_i^a,(nu+alpha_i^2)/(nu-1)*\tilde\Sigma_i),
# %   qb(i) = c_{beta_i}/sigma_i*L_{n-1,nu-1}(a(i),b(i),\tilde\mu_i^b,(nu+beta_i^2)/(nu-1)*\tilde\Sigma_i),
# %   where c_x = gamma((nu-1)/2)*sqrt(nu)/(2*sqrt(pi)*gamma(nu/2))*(1+x^2/nu)^(-(nu-1)/2).
# %   The program also deals with the cases with nu<=1.
# %

qfunT = function(a,b,S,nu){
  # print("qfunT")
  # print(nu)
  n = length(a)
  if(nu<=-1){
    qa = qb = matrix(NaN,n,1)
    return(list(qa = qa,qb = qb))
  }
  s = sqrt(diag(as.matrix(S)))
  as = a/s
  bs = b/s
  if(n==1){
    if(nu!=1){
      cc = gamma((nu-1)/2)*sqrt(nu)/(2*sqrt(pi)*gamma(nu/2))
      qa = cc/s/(1+as^2/nu)^((nu-1)/2)
      qb = cc/s/(1+bs^2/nu)^((nu-1)/2)
    }else{
      if(nu>0){
        qa = -log(1+as^2)/(2*pi*s)
        qb = -log(1+bs^2)/(2*pi*s)
      }else{
        qa = NaN
        qb = NaN
      }
    }
    if(is.infinite(qa)){
      qa = 0
    }
    if(is.infinite(qb)){
      qb = 0
    }
    return(list(qa = qa,qb = qb))
  }
  qa = qb = matrix(0,n,1)

  for(i in 1:n){
    RR = S[-i,-i] - S[-i,i]%*%t(S[-i,i])/S[i,i]
    ma = S[-i,i]/S[i,i]*a[i]
    a1 = a[-i]-ma
    b1 = b[-i]-ma
    mb = S[-i,i]/S[i,i]*b[i]
    a2 = a[-i]-mb
    b2 = b[-i]-mb
    if(nu==1){
      m = sum(is.infinite(a[-i])|is.infinite(b[-i]))
    }

    if(is.finite(a[i])){
      da = nu+as[i]^2
      if(nu==1 && m==n-1){
        #% This is the case that qa(i) and qb(i) do not exist but qa(i)-qb(i) exists
        if(is.finite(b[i])){
          db = nu+bs[i]^2
          #% For n=2, we have explicit expression from note05
          if(n==2){

            qa[i] = log(db/da)
            if(is.finite(a1)){
              qa[i] = qa[i]/2-asinh(a1/sqrt(RR)/sqrt(da))+asinh(a2/sqrt(RR)/sqrt(db))
            }
            if(is.finite(b1)){
              qa[i] = qa[i]/2+asinh(b1/sqrt(RR)/sqrt(da))-asinh(b2/sqrt(RR)/sqrt(db))
            }

          }else{
            # if(n<=4){
            #   f1 = function(t,z,c1,d1){
            #     #%  Special case for nu=1, which requires direct integration
            #     st = sqrt(t)
            #     return(pmvnorm(lower = st*c1,upper = st*d1,sigma = RR)*exp(-z*t/2)/t)
            #   }
            #   f12 = function(t){
            #     return(sapply(X = t,FUN = f1,z = da,c1 = a1,d1 = b1)-sapply(X = t,FUN = f1,z = db,c1 = a2,d1 = b2))
            #   }
            #   qa[i] = integrate(f12,0,Inf,abs.tol = 1e-10,rel.tol = 1e-10)$value
            # }else{
            nu1 = nu+1e-1
            nu2 = nu-1e-1
            # nu1 = nu+1
            # nu2 = nu-1
            qa1 = lfun(a1,b1,RR*(nu1+as[i]^2),nu1-1)/(1+as[i]^2/nu1)^((nu1-1)/2)-
              lfun(a2,b2,RR*(nu1+bs[i]^2),nu1-1)/(1+bs[i]^2/nu1)^((nu1-1)/2)
            qa2 = lfun(a1,b1,RR*(nu2+as[i]^2),nu2-1)/(1+as[i]^2/nu2)^((nu2-1)/2)-
              lfun(a2,b2,RR*(nu2+bs[i]^2),nu2-1)/(1+bs[i]^2/nu2)^((nu2-1)/2)
            qa[i] = (qa1+qa2)/2
            #}
          }
        }
      }else{
        qa[i] = lfun(a1,b1,RR*da,nu-1)/(da/nu)^((nu-1)/2)
      }
    }

    if(is.finite(b[i])){
      db = nu+bs[i]^2
      if(nu!=1||m!=n-1){
        qb[i] = lfun(a2,b2,RR*db,nu-1)/(db/nu)^((nu-1)/2)
      }
    }
  }
  cc = sqrt(nu)/(2*sqrt(pi)*gamma(nu/2))
  qa = cc*qa/s
  qb = cc*qb/s
  return(list(qa = qa, qb = qb))
}


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

qfunT_a = function(a,mu,S,nu){
  n = length(mu)
  s = sqrt(diag(as.matrix(S)))
  as = (a-mu)/s
  q = sqrt(nu)*gamma((nu-1)/2)/(2*sqrt(pi)*gamma(nu/2))*(1+as^2/nu)^(-(nu-1)/2)/s
  if(n>1){
    for(i in 1:n){
      if(is.finite(a[i])){
        muj  = mu[-i]+S[-i,i]*as[i]/s[i]
        Sj   = (nu+as[i]^2)/(nu-1)*(S[-i,-i]-S[-i,i]%*%t(S[i,-i])/S[i,i])
        sj   = sqrt(diag(as.matrix(Sj)))
        Rj   = Sj/(sj%*%t(sj))
        q[i] = log2prod(q[i],pmvt.genz(upper = -a[-i],mean = -muj,nu = nu-1,sigma = Sj,uselog2 = TRUE)$Estimation)
      }
    }
  }
  return(q)
}
