# #########################################################################################################
# #########################################################################################################
# 
withinfsT = function(lower = rep(-Inf,length(mu)),upper = rep(Inf,length(mu)),mu,Sigma,nu){
  bool = is.infinite(lower) & is.infinite(upper)
  p          = length(bool)
  seqq       = seq(1,p)
  bool       = seqq[bool]
  q          = length(bool) #infinite dims
  
  
  op2 = list()
  op2$mean = matrix(NA,p,1)
  
  out.finite = meanvarTall(lower[-bool],upper[-bool],mu[-bool],Sigma[-bool,-bool],nu,omega = TRUE)
  
  op2$mean[bool]  = c(mu[bool] + Sigma[bool,-bool]%*%solve(Sigma[-bool,-bool])%*%(out.finite$mean - mu[-bool]))
  op2$mean[-bool] = out.finite$mean
  op2$EYY    = matrix(NA,p,p)
  op2$varcov = matrix(NA,p,p)
  p2 = p - q
  
  #For nu>2
  #omega12 = nu/(nu-2)*pmvnormt(lower = lower,upper = upper,mean = mu,sigma = nu*Sigma/(nu-2),nu = nu - 2)/pmvnormt(lower = lower,upper = upper,mean = mu,sigma = Sigma,nu = nu)
  
  op2$varcov[bool,bool]  = out.finite$omega*Sigma[bool,bool] - 
                           Sigma[bool,-bool]%*%
                            solve(Sigma[-bool,-bool])%*%
                            (out.finite$omega*diag(p2) - 
                               out.finite$varcov%*%solve(Sigma[-bool,-bool]))%*%Sigma[-bool,bool]
  
  op2$varcov[-bool,-bool] = out.finite$varcov
  op2$varcov[bool,-bool] = Sigma[bool,-bool]%*%solve(Sigma[-bool,-bool])%*%out.finite$varcov
  op2$varcov[-bool,bool] = t(op2$varcov[bool,-bool])
  
  #-----------------------------------------------
  ind = seq_len(p)[is.infinite(lower)|is.infinite(upper)]
  d = p - length(ind) + nu
  
  #Check existence of moments
  if(nu<=2){
    if(d<=1){
      op2$varcov[ind,] = NaN
      op2$varcov[,ind] = NaN
      op2$mean[ind]   = NaN
    }else{
      if(d<=2){
        op2$varcov[ind,ind] = NaN
      }
    }
  }
  #-----------------------------------------------
  
  op2$EYY = op2$varcov + op2$mean%*%t(op2$mean)
  
  return(op2)
}
# 
# 
# #########################################################################################################
# 

withinfs_meanT = function(lower = rep(-Inf,length(mu)),upper = rep(Inf,length(mu)),mu,Sigma,nu){
  bool = is.infinite(lower) & is.infinite(upper)
  p = length(bool)
  seqq   = seq(1,p)
  bool   = seqq[bool]
  q      = length(bool) #infinite dims
  out.finite = onlymeanT(lower[-bool],upper[-bool],mu[-bool],Sigma[-bool,-bool],nu)
  
  op2 = list()
  op2$mean = matrix(NA,p,1)
  op2$mean[bool]  = c(mu[bool] + Sigma[bool,-bool]%*%solve(Sigma[-bool,-bool])%*%(out.finite$mean - mu[-bool]))
  op2$mean[-bool] = out.finite$mean
  
  #-----------------------------------------------
  ind = seq_len(p)[is.infinite(lower)|is.infinite(upper)]
  d = p - length(ind) + nu
  
  #Check existence of moments
  if(nu<=2 & d<=1){
      op2$mean[ind]   = NaN
  }
  #-----------------------------------------------
  
  return(op2)
}

#########################################################################################################

withinfs_meanT0 = function(lower = rep(-Inf,length(mu)),upper = rep(Inf,length(mu)),mu,Sigma,nu){
  bool = is.infinite(lower) & is.infinite(upper)
  p = length(bool)
  seqq   = seq(1,p)
  bool   = seqq[bool]
  q      = length(bool) #infinite dims
  out.finite = onlymeanT0(lower[-bool],upper[-bool],mu[-bool],Sigma[-bool,-bool],nu)
  
  op2 = list()
  op2$mean = matrix(NA,p,1)
  op2$mean[bool]  = c(mu[bool] + Sigma[bool,-bool]%*%solve(Sigma[-bool,-bool])%*%(out.finite$mean - mu[-bool]))
  op2$mean[-bool] = out.finite$mean
  
  #-----------------------------------------------
  ind = seq_len(p)[is.infinite(lower)|is.infinite(upper)]
  d = p - length(ind) + nu
  
  #Check existence of moments
  if(nu<=2 & d<=1){
    op2$mean[ind]   = NaN
  }
  #-----------------------------------------------
  
  return(list(mean = op2,logF00 = out.finite$logF0))
}

#########################################################################################################

# 
# library(tlrmvnmvt)
# 
# MCT = function(n,a,b,mu,S,nu,algo = "rejection"){
#   gen = tmvtnorm::rtmvt(n = n,mean = mu,sigma = S,df = nu,lower = a,upper = b,algorithm = algo)
#   return(list(mean = colMeans(gen),varcov = var(gen)))
# }
# 
# #########################################################################################################
# nu = 4
# #########################################################################################################
# 
# set.seed(0)
# p = 4
# mu  = c(1:p/10)
# lower = a = -mu - abs(rnorm(p))
# upper = b = mu+1 +0.5*a + abs(rnorm(p)) + 3
# s  = matrix(rnorm(p^2),p,p)
# S = round(s%*%t(s),2)
# Sigma = S
# 
# 
# lower[(1:p)<(p/2)] = -Inf
# upper[(1:p)<(p/2)] = Inf
# 
# 
# bool = is.infinite(lower) & is.infinite(upper)
# 
# 
# meanvarT.Lin.LRIC(lower,upper,mu,Sigma,nu,TRUE)
# meanvarTMD(lower,upper,mu,Sigma,nu = nu,dist = "t")
# 
# # mean.cond = withinfs_meanT(lower,upper,mu,Sigma,nu,bool)$mean
# # mean.all  = onlymeanT(lower,upper,mu,Sigma,nu)$mean
# # mean.MC = MCT(n = 200000,lower,upper,mu,Sigma,nu)$mean
# # 
# # cbind(mean.cond,mean.all,mean.MC)
# 
# meanvar.cond = withinfsTold(lower,upper,mu,Sigma,nu,bool)
# meanvar.all  = meanvarT16(lower,upper,mu,Sigma,nu)
# meanvar.MC = MCT(n = 100000,lower,upper,mu,Sigma,nu)
# 
# cbind(meanvar.cond$mean,meanvar.all$mean,meanvar.MC$mean)
# meanvar.cond$varcov;meanvar.all$varcov;meanvar.MC$varcov
# 
# ###########################################
# 
# lower = a = -mu - abs(rnorm(p))
# upper = b = mu+1 +0.5*a + abs(rnorm(p)) + 3
# 
# lower[(1:p)<(p/2 +1)] = -Inf
# upper[(1:p)<(p/2 +1)] = Inf
# 
# bool = is.infinite(lower) & is.infinite(upper)
# 
# # mean.cond = withinfs_meanT(lower,upper,mu,Sigma,nu,bool)$mean
# # mean.all  = onlymeanT(lower,upper,mu,Sigma,nu)$mean
# # mean.MC = MCT(n = 200000,lower,upper,mu,Sigma,nu)$mean
# # 
# # cbind(mean.cond,mean.all,mean.MC)
# 
# meanvar.cond = withinfsTold(lower,upper,mu,Sigma,nu,bool)
# meanvar.all  = meanvarT16(lower,upper,mu,Sigma,nu)
# meanvar.MC = MCT(n = 100000,lower,upper,mu,Sigma,nu)
# 
# cbind(meanvar.cond$mean,meanvar.all$mean,meanvar.MC$mean)
# meanvar.cond$varcov;meanvar.all$varcov;meanvar.MC$varcov
# 
# 
# ###########################################
# 
# lower = a = -mu - abs(rnorm(p))
# upper = b = mu+1 +0.5*a + abs(rnorm(p)) + 3
# 
# lower[(1:p)<(p/2 +2)] = -Inf
# upper[(1:p)<(p/2 +2)] = Inf
# 
# bool = is.infinite(lower) & is.infinite(upper)
# 
# # mean.cond = withinfs_meanT(lower,upper,mu,Sigma,nu,bool)$mean
# # mean.all  = onlymeanT(lower,upper,mu,Sigma,nu)$mean
# # mean.MC = MCT(n = 200000,lower,upper,mu,Sigma,nu)$mean
# # 
# # cbind(mean.cond,mean.all,mean.MC)
# 
# meanvar.cond = withinfsTold(lower,upper,mu,Sigma,nu,bool)
# meanvar.all  = meanvarT16(lower,upper,mu,Sigma,nu)
# meanvar.MC = MCT(n = 100000,lower,upper,mu,Sigma,nu)
# 
# cbind(meanvar.cond$mean,meanvar.all$mean,meanvar.MC$mean)
# meanvar.cond$varcov;meanvar.all$varcov;meanvar.MC$varcov