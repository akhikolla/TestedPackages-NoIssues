RcppMCT.lin = function(n,a,b,mu,S,nu = 10000,omega = FALSE){
  gen = RcppTT.GS(n=n, mu, as.matrix(S), nu, lower=a, upper=b)
  mean0 = as.matrix(colMeans(gen))
  varcov = var(gen)
  
  #-----------------------------------------------
  p = length(mu)
  ind = seq_len(p)[is.infinite(a)|is.infinite(b)]
  d = p - length(ind) + nu
  
  #Check existence of moments
  if(nu<=2){
    if(d<=1){
      varcov[ind,] = NaN
      varcov[,ind] = NaN
      mean0[ind]   = NaN
    }else{
      if(d<=2){
        varcov[ind,ind] = NaN
      }
    }
  }
  #-----------------------------------------------
  
  if(omega){
    dd = mahalanobis(x = gen,center = mu,cov = as.matrix(S))
    omega12 = (nu + mean(dd))/(nu + p - 2)
    return(list(mean = mean0,EYY = varcov + mean0%*%t(mean0), varcov = varcov,omega = omega12))
  }else{
    return(list(mean = mean0,EYY = varcov + mean0%*%t(mean0), varcov = varcov))
  }
}


RcppTT.GS = function(n, mu=rep(0,nrow(S)), S=diag(length(mu)), nu=2, lower=rep(-Inf, length(mu)), upper=rep(Inf, length(mu)))
  
{
  
  #  require(mvtnorm)
  
  p=length(mu)
  
  ## Verify error at parameters specification
  # if(length(lower) != length(upper)) stop("The lengths of the 'lower' and 'upper' truncated values must be equal!")
  # for(i in 1:length(lower))
  # {
  #   if(upper[i] <= lower[i]) stop("The lower limit is larger than upper limit for truncation!")
  # }
  # if(length(lower) != nrow(S) | length(upper) != nrow(S) ) stop("The dimension of scale-covariance matrix must be equal to the length of the 'lower/upper' truncated value!")
  # if(nu <=0) stop("The degree of freedom must be larger than zero!")
  # if(det(S)<=0) stop("The S matrix must be inversible!")
  # 
  
  s = sqrt(diag(S))
  
  R = S/outer(s,s,"*")
  
  x=qt(runif(rep(1,p),pt((lower-mu)/s,df=nu),pt((upper-mu)/s,df=nu)),df=nu)
  a=ifelse(lower==-Inf,rep(-1e12,p),(lower-mu)/s)
  b=ifelse(lower== Inf,rep( 1e12,p),(upper-mu)/s)
  Z = TT_GS_sp(n,R,nu,x,lower=a,upper=b)
  
  X = t(mu + t(Z) * s)
  
  return(X)
  
}

# MCT = function(n,a,b,mu,S,nu,algo = "rejection",omega = FALSE){
#   gen = tmvtnorm::rtmvt(n = n,mean = mu,sigma = as.matrix(S),df = nu,lower = a,upper = b,algorithm = algo)
#   mean0 = as.matrix(colMeans(gen))
#   varcov = var(gen)
#   
#   #-----------------------------------------------
#   p = length(mu)
#   ind = seq_len(p)[is.infinite(a)|is.infinite(b)]
#   d = p - length(ind) + nu
#   
#   #Check existence of moments
#   if(nu<=2){
#     if(d<=1){
#       varcov[ind,] = NaN
#       varcov[,ind] = NaN
#       mean0[ind]   = NaN
#     }else{
#       if(d<=2){
#         varcov[ind,ind] = NaN
#       }
#     }
#   }
#   #-----------------------------------------------
#   
#   if(omega){
#   dd = mahalanobis(x = gen,center = mu,cov = as.matrix(S))
#   omega12 = (nu + mean(dd))/(nu + p - 2)
#   return(list(mean = mean0,EYY = varcov + mean0%*%t(mean0), varcov = varcov,omega = omega12))
#   }else{
#     return(list(mean = mean0,EYY = varcov + mean0%*%t(mean0), varcov = varcov))
#   }
#   
# }

# MCT.lin = function(n,a,b,mu,S,nu,omega = FALSE){
#   gen = TT.GS(n=n, mu, as.matrix(S), nu, lower=a, upper=b)
#   mean0 = as.matrix(colMeans(gen))
#   varcov = var(gen)
#   
#   #-----------------------------------------------
#   p = length(mu)
#   ind = seq_len(p)[is.infinite(a)|is.infinite(b)]
#   d = p - length(ind) + nu
#   
#   #Check existence of moments
#   if(nu<=2){
#     if(d<=1){
#       varcov[ind,] = NaN
#       varcov[,ind] = NaN
#       mean0[ind]   = NaN
#     }else{
#       if(d<=2){
#         varcov[ind,ind] = NaN
#       }
#     }
#   }
#   #-----------------------------------------------
#   
#   if(omega){
#     dd = mahalanobis(x = gen,center = mu,cov = as.matrix(S))
#     omega12 = (nu + mean(dd))/(nu + p - 2)
#     return(list(mean = mean0,EYY = varcov + mean0%*%t(mean0), varcov = varcov,omega = omega12))
#   }else{
#     return(list(mean = mean0,EYY = varcov + mean0%*%t(mean0), varcov = varcov))
#   }
# }

#


#Slice sampling from momentTT package

# TT.GS = function(n, mu=rep(0,nrow(S)), S=diag(length(mu)), nu=2, lower=rep(-Inf, length(mu)), upper=rep(Inf, length(mu)))
#   
# {
#   
#   #  require(mvtnorm)
#   
#   p=length(mu)
#   
#   ## Verify error at parameters specification
#   if(length(lower) != length(upper)) stop("The lengths of the 'lower' and 'upper' truncated values must be equal!")
#   for(i in 1:length(lower))
#   {
#     if(upper[i] <= lower[i]) stop("The lower limit is larger than upper limit for truncation!")
#   }
#   if(length(lower) != nrow(S) | length(upper) != nrow(S) ) stop("The dimension of scale-covariance matrix must be equal to the length of the 'lower/upper' truncated value!")
#   if(nu <=0) stop("The degree of freedom must be larger than zero!")
#   if(det(S)<=0) stop("The S matrix must be inversible!")
#   
#   
#   s = sqrt(diag(S))
#   
#   R = S/outer(s,s,"*")
#   
#   Z = TT.GS.sp(n,R,nu,lower=(lower-mu)/s,upper=(upper-mu)/s)
#   
#   X = t(mu + t(Z) * s)
#   
#   return(X)
#   
# }

# TT.GS.sp = function(n,R,nu,lower,upper)
#   
# {
#   
#   #initial value
#   x=qt(runif(rep(1,p),pt(lower,df=nu),pt(upper,df=nu)),df=nu)
#   
#   if(n<1) return(t(x))    
#   
#   X = matrix(NA, n, p)
#   
#   R.inv = solve(R)
#   
#   for(i in 1:n)
#     
#   {
#     
#     delta = sum(colSums(x*R.inv)*x)
#     
#     y = runif(1,0,exp(-.5*(nu+p)*log(1+delta/nu)))
#     
#     kap = nu*(y^(-2/(nu+p))-1)
#     
#     for(j in 1:p)
#       
#     {
#       
#       ss = x[-j]%*%R.inv[-j,-j]%*%x[-j]
#       
#       mj = - sum(R.inv[-j,j]*x[-j]) / R.inv[j,j]
#       
#       tj = sqrt(mj^2 + (kap - ss) / R.inv[j,j])
#       
#       xij = runif(1,max(lower[j],mj-tj),min(upper[j],mj+tj))
#       
#       X[i,j] = xij
#       
#       x[j] = xij
#       
#     }
#     
#   }
#   
#   return(X)
#   
# }