#########################################################################################################
#########################################################################################################

corrector = function(lower = rep(-Inf,length(mu)),upper = rep(Inf,length(mu)),mu,Sigma,bw=36){
  p = length(lower)
  ss = sqrt(diag(as.matrix(Sigma)))
  liminf = mu - bw*ss
  limsup = mu + bw*ss
  bool1 = upper < liminf
  bool2 = lower > limsup
  while(sum(bool1) +  sum(bool2) == 0){
    bw = bw - 2
    liminf = mu - bw*ss
    limsup = mu + bw*ss
    bool1 = upper < liminf
    bool2 = lower > limsup
  }
  bool   = bool1 | bool2
  seqq   = seq(1,p)
  bool   = seqq[bool]
  
  bool1  = seqq[bool1]
  bool2  = seqq[bool2]
  
  val    = rep(NA,p)
  vars   = rep(NA,p)
  for(i in bool){
    #print(i)
    mv1 = meanvarNuni(lower[i],upper[i],mu[i],ss[i])
    val[i] = mv1$mean
    vars[i] = mv1$varcov
  }
  val  = val[bool]
  vars = vars[bool]
  
  if(length(bool) == p){
    out1 = matrix(val,p,1)
    out3 = matrix(1e-10,p,p)
    diag(out3) = vars
    out2 = out3 + out1%*%t(out1)
    return(list(mean = out1,EYY = out2,varcov = out3))
  }
  tilSj = Sigma[-bool,-bool] - Sigma[-bool,bool]%*%solve(Sigma[bool,bool])%*%Sigma[bool,-bool]
  op1 = meanvarN7(lower = lower[-bool],
                  upper = upper[-bool],
                  mu = c(mu[-bool] + Sigma[-bool,bool]%*%solve(Sigma[bool,bool])%*%(val - mu[bool])),
                  Sigma = tilSj)
  op2 = op1
  op2$mean = matrix(NA,p,1)
  op2$mean[bool]  = val
  op2$mean[-bool] = op1$mean
  op2$varcov = matrix(1e-10,p,p)
  op2$varcov[bool,bool] = vars
  op2$varcov[-bool,-bool] = op1$varcov
  op2$EYY = op2$varcov + op2$mean%*%t(op2$mean)
  return(op2)
}

#corrector(upper = c(0,-50,-50),mu=mu,Sigma=Sigma)

#########################################################################################################

corrector_onlymean = function(lower = rep(-Inf,length(mu)),upper = rep(Inf,length(mu)),mu,Sigma,bw=36){
  p = length(lower)
  ss = sqrt(diag(as.matrix(Sigma)))
  liminf = mu - bw*ss
  limsup = mu + bw*ss
  bool1 = upper < liminf
  bool2 = lower > limsup
  
  while(sum(bool1) +  sum(bool2) == 0){
    bw = bw - 2
    liminf = mu - bw*ss
    limsup = mu + bw*ss
    bool1 = upper < liminf
    bool2 = lower > limsup
  }
  bool   = bool1 | bool2
  seqq   = seq(1,p)
  bool   = seqq[bool]
  
  bool1  = seqq[bool1]
  bool2  = seqq[bool2]
  
  val    = rep(NA,p)
  for(i in bool){
    #print(i)
    mv1 = onlymeanNuni(lower[i],upper[i],mu[i],ss[i])
    val[i] = mv1$mean
  }
  val  = val[bool]
  
  if(length(bool) == p){
    out1 = matrix(val,p,1)
    return(list(mean = out1))
  }
  tilSj = Sigma[-bool,-bool] - Sigma[-bool,bool]%*%solve(Sigma[bool,bool])%*%Sigma[bool,-bool]
  op1 = Vaida.LRIC.onlymean(a = lower[-bool],
                            b = upper[-bool],
                            mu = c(mu[-bool] + Sigma[-bool,bool]%*%solve(Sigma[bool,bool])%*%(val - mu[bool])),
                            Sigma = tilSj)
  op2 = op1
  op2$mean = matrix(NA,p,1)
  op2$mean[bool]  = val
  op2$mean[-bool] = op1$mean
  return(op2)
}

#########################################################################################################

withinfs = function(lower = rep(-Inf,length(mu)),upper = rep(Inf,length(mu)),mu,Sigma,bool){
  #bool = is.infinite(lower) & is.infinite(upper)
  p = length(bool)
  seqq   = seq(1,p)
  bool   = seqq[bool]
  q      = length(bool) #infinite dims
  if(q < 10){
    out.finite = Kan.LRIC(lower[-bool],upper[-bool],mu[-bool],Sigma[-bool,-bool])
  }else{
    out.finite = Vaida.LRIC(lower[-bool],upper[-bool],mu[-bool],Sigma[-bool,-bool])
  }
  op2 = list()
  op2$mean = matrix(NA,p,1)
  op2$mean[bool]  = c(mu[bool] + Sigma[bool,-bool]%*%solve(Sigma[-bool,-bool])%*%(out.finite$mean - mu[-bool]))
  op2$mean[-bool] = out.finite$mean
  op2$EYY    = matrix(NA,p,p)
  op2$varcov = matrix(NA,p,p)
  op2$varcov[bool,bool]  = Sigma[bool,bool] - Sigma[bool,-bool]%*%solve(Sigma[-bool,-bool])%*%(diag(p-q) - out.finite$varcov%*%solve(Sigma[-bool,-bool]))%*%Sigma[-bool,bool]
  op2$varcov[-bool,-bool] = out.finite$varcov
  op2$varcov[bool,-bool] = Sigma[bool,-bool]%*%solve(Sigma[-bool,-bool])%*%out.finite$varcov
  op2$varcov[-bool,bool] = t(op2$varcov[bool,-bool])
  op2$EYY = op2$varcov + op2$mean%*%t(op2$mean)
  return(op2)
}

#########################################################################################################

withinfs_mean = function(lower = rep(-Inf,length(mu)),upper = rep(Inf,length(mu)),mu,Sigma,bool){
  #bool = is.infinite(lower) & is.infinite(upper)
  p = length(bool)
  seqq   = seq(1,p)
  bool   = seqq[bool]
  q      = length(bool) #infinite dims
  out.finite = Vaida.LRIC.onlymean(lower[-bool],upper[-bool],mu[-bool],Sigma[-bool,-bool])
  
  op2 = list()
  op2$mean = matrix(NA,p,1)
  op2$mean[bool]  = c(mu[bool] + Sigma[bool,-bool]%*%solve(Sigma[-bool,-bool])%*%(out.finite$mean - mu[-bool]))
  op2$mean[-bool] = out.finite$mean
  return(op2)
}
