# #######################################################################
# #STUDENT T
# #######################################################################
# 
# meanvarT = function(a,b,mu,Sigma,nu)
# {
#   #if(nu>=4){
#   p = length(mu)
#   nnu = nu/(nu-2)
#   if(p==1){
#     F0 = pent(b,mu,Sigma,nu) - pent(a,mu,Sigma,nu)
#     nnusigma2 = nnu*Sigma
#     ta = dent(a,mu,nnusigma2,nu-2)
#     tb = dent(b,mu,nnusigma2,nu-2)
#     F1 = mu*F0 + nnusigma2*(ta-tb)
#     F2 = mu*F1 + nnusigma2*(pent(b,mu,nnusigma2,nu-2) - pent(a,mu,nnusigma2,nu-2) + ifelse(a==-Inf,0,a*ta) - ifelse(b==Inf,0,b*tb))
#     return(list(mean = F1/F0,EYY = F2/F0,varcov = F2/F0 - (F1/F0)^2))
#   }
#   #GB = GenzBretz(maxpts = (p-1)*1e4, abseps = 1e-6, releps = 0)
#   #print(GB$maxpts)
#   #F0 = pmvt(lower = a-mu,upper = b-mu,df = nu,sigma = Sigma)[1]
# 
#   logF0 = pmvt.genz(lower = a-mu,upper = b-mu,nu = nu,sigma = Sigma,uselog2 = TRUE,N = 799)$Estimation
#   F0nnu = pmvt.genz(lower = a-mu,upper = b-mu,nu = nu - 2,sigma = nnu*Sigma,N = 799)$Estimation
# 
#   #Vectors ca and cb
#   SSigma  = nnu*Sigma
#   ssigma2 = diag(SSigma)
#   ca = cb = a0 = a1 = rep(0,p)
#   deltaA  = (nu - 2 + ((a - mu)^2)/diag(SSigma))/(nu - 1)
#   deltaB  = (nu - 2 + ((b - mu)^2)/diag(SSigma))/(nu - 1)
#   yA      = (a - mu)/diag(SSigma)
#   yB      = (b - mu)/diag(SSigma)
# 
#   Wa = Wb = matrix(0,p,p)
# 
#   for(j in 1:p)
#   {
#     #W matrix construction
# 
#     if(a[j]!=-Inf){
#       aux1     = onlymeanT0(a = a[-j],b = b[-j],mu = mu[-j] + yA[j]*SSigma[,j][-j],Sigma = deltaA[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),nu = nu-1)
#       ca[j]    = log2prod(dent(a[j],mu[j],ssigma2[j],nu-2),aux1$logF00)
#       Wa[-j,j] = aux1$mean
#       Wa[j,j]  = a[j]
#     }
#     if(b[j]!= Inf){
#       aux2     = onlymeanT0(a = a[-j],b = b[-j],mu = mu[-j] + yB[j]*SSigma[,j][-j],Sigma = deltaB[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),nu = nu-1)
#       cb[j]    = log2prod(dent(b[j],mu[j],ssigma2[j],nu-2),aux2$logF00)
#       Wb[-j,j] = aux2$mean
#       Wb[j,j]  = b[j]
#     }
#   }
#   muY  = mu + log2ratio(SSigma%*%(ca - cb),logF0)
#   Exx  = muY%*%t(mu) +  log2ratio((F0nnu*diag(p) + Wa%*%diag(ca) - Wb%*%diag(cb))%*%SSigma,logF0)
#   Exx = (Exx + t(Exx))/2
#   varY = Exx - muY%*%t(muY)
#   return(list(mean = muY,EYY = Exx,varcov = varY))
# }
# 
# #######################################################################
# #######################################################################
# #######################################################################
# 
# meanvarT_lower = function(a,mu,Sigma,nu)
# {
#   #if(nu>=4){
#   p = length(mu)
#   nnu = nu/(nu-2)
#   if(p==1){
#     F0 = 1 - pent(a,mu,Sigma,nu)
#     nnusigma2 = nnu*Sigma
#     ta = dent(a,mu,nnusigma2,nu-2)
#     tb = 0
#     F1 = mu*F0 + nnusigma2*(ta-tb)
#     F2 = mu*F1 + nnusigma2*(1 - pent(a,mu,nnusigma2,nu-2) + ifelse(a==-Inf,0,a*ta))
#     return(list(mean = F1/F0,EYY = F2/F0,varcov = F2/F0 - (F1/F0)^2))
#   }
#   #GB = GenzBretz(maxpts = (p-1)*1e4, abseps = 1e-6, releps = 0)
#   #print(GB$maxpts)
#   #F0 = pmvt(lower = a-mu,upper = b-mu,df = nu,sigma = Sigma)[1]
# 
#   logF0 = pmvt.genz(lower = a-mu,nu = nu,sigma = Sigma,uselog2 = TRUE,N = 799)$Estimation
#   F0nnu = pmvt.genz(lower = a-mu,nu = nu - 2,sigma = nnu*Sigma,N = 799)$Estimation
# 
#   #Vectors ca and cb
#   SSigma  = nnu*Sigma
#   ssigma2 = diag(SSigma)
#   ca = a0 = a1 = rep(0,p)
#   deltaA  = (nu - 2 + ((a - mu)^2)/diag(SSigma))/(nu - 1)
#   #deltaB  = (nu - 2 + ((b - mu)^2)/diag(SSigma))/(nu - 1)
#   yA      = (a - mu)/diag(SSigma)
#   #yB      = (b - mu)/diag(SSigma)
# 
#   Wa = matrix(0,p,p)
# 
#   for(j in 1:p)
#   {
#     #W matrix construction
# 
#     if(a[j]!=-Inf){
#       aux1     = onlymeanT0(a = a[-j],b = rep(Inf,p-1),mu = mu[-j] + yA[j]*SSigma[,j][-j],Sigma = deltaA[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),nu = nu-1)
#       ca[j]    = log2prod(dent(a[j],mu[j],ssigma2[j],nu-2),aux1$logF00)
#       Wa[-j,j] = aux1$mean
#       Wa[j,j]  = a[j]
#     }
#     # if(b[j]!= Inf){
#     #   aux2     = onlymeanT0(a = a[-j],b = b[-j],mu = mu[-j] + yB[j]*SSigma[,j][-j],Sigma = deltaB[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),nu = nu-1)
#     #   cb[j]    = log2prod(dent(b[j],mu[j],ssigma2[j],nu-2),aux2$logF00)
#     #   Wb[-j,j] = aux2$mean
#     #   Wb[j,j]  = b[j]
#     # }
#   }
#   muY  = mu + log2ratio(SSigma%*%ca,logF0)
#   Exx  = muY%*%t(mu) +  log2ratio((F0nnu*diag(p) + Wa%*%diag(ca))%*%SSigma,logF0)
#   Exx = (Exx + t(Exx))/2
#   varY = Exx - muY%*%t(muY)
#   return(list(mean = muY,EYY = Exx,varcov = varY))
# }
# 
# #######################################################################
# #######################################################################
# #######################################################################
# 
# meanvarT_upper = function(b,mu,Sigma,nu)
# {
#   #if(nu>=4){
#   p = length(mu)
#   nnu = nu/(nu-2)
#   if(p==1){
#     F0 = pent(b,mu,Sigma,nu)
#     nnusigma2 = nnu*Sigma
#     tb = dent(b,mu,nnusigma2,nu-2)
#     F1 = mu*F0 - nnusigma2*tb
#     F2 = mu*F1 + nnusigma2*(pent(b,mu,nnusigma2,nu-2) - ifelse(b==Inf,0,b*tb))
#     return(list(mean = F1/F0,EYY = F2/F0,varcov = F2/F0 - (F1/F0)^2))
#   }
#   #GB = GenzBretz(maxpts = (p-1)*1e4, abseps = 1e-6, releps = 0)
#   #print(GB$maxpts)
#   #F0 = pmvt(lower = a-mu,upper = b-mu,df = nu,sigma = Sigma)[1]
# 
#   logF0 = pmvt.genz(upper = b-mu,nu = nu,sigma = Sigma,uselog2 = TRUE,N = 799)$Estimation
#   F0nnu = pmvt.genz(upper = b-mu,nu = nu - 2,sigma = nnu*Sigma,N = 799)$Estimation
# 
#   #Vectors ca and cb
#   SSigma  = nnu*Sigma
#   ssigma2 = diag(SSigma)
#   cb = a0 = a1 = rep(0,p)
#   #deltaA  = (nu - 2 + ((a - mu)^2)/diag(SSigma))/(nu - 1)
#   deltaB  = (nu - 2 + ((b - mu)^2)/diag(SSigma))/(nu - 1)
#   #yA      = (a - mu)/diag(SSigma)
#   yB      = (b - mu)/diag(SSigma)
# 
#   Wb = matrix(0,p,p)
# 
#   for(j in 1:p)
#   {
#     #W matrix construction
# 
#     # if(a[j]!=-Inf){
#     #   aux1     = onlymeanT0(a = a[-j],b = b[-j],mu = mu[-j] + yA[j]*SSigma[,j][-j],Sigma = deltaA[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),nu = nu-1)
#     #   ca[j]    = log2prod(dent(a[j],mu[j],ssigma2[j],nu-2),aux1$logF00)
#     #   Wa[-j,j] = aux1$mean
#     #   Wa[j,j]  = a[j]
#     # }
#     if(b[j]!= Inf){
#       aux2     = onlymeanT0(a =rep(-Inf,p-1),b = b[-j],mu = mu[-j] + yB[j]*SSigma[,j][-j],Sigma = deltaB[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),nu = nu-1)
#       cb[j]    = log2prod(dent(b[j],mu[j],ssigma2[j],nu-2),aux2$logF00)
#       Wb[-j,j] = aux2$mean
#       Wb[j,j]  = b[j]
#     }
#   }
#   muY  = mu - log2ratio(SSigma%*%cb,logF0)
#   Exx  = muY%*%t(mu) +  log2ratio((F0nnu*diag(p) - Wb%*%diag(cb))%*%SSigma,logF0)
#   Exx = (Exx + t(Exx))/2
#   varY = Exx - muY%*%t(muY)
#   return(list(mean = muY,EYY = Exx,varcov = varY))
# }
# 
# 
# #######################################################################
# #######################################################################
# #######################################################################
# 
# meanvarT_finite = function(a,b,mu,Sigma,nu)
# {
#   #if(nu>=4){
#   p = length(mu)
#   nnu = nu/(nu-2)
#   if(p==1){
#     F0 = pent(b,mu,Sigma,nu) - pent(a,mu,Sigma,nu)
#     nnusigma2 = nnu*Sigma
#     ta = dent(a,mu,nnusigma2,nu-2)
#     tb = dent(b,mu,nnusigma2,nu-2)
#     F1 = mu*F0 + nnusigma2*(ta-tb)
#     F2 = mu*F1 + nnusigma2*(pent(b,mu,nnusigma2,nu-2) - pent(a,mu,nnusigma2,nu-2) + a*ta - b*tb)
#     return(list(mean = F1/F0,EYY = F2/F0,varcov = F2/F0 - (F1/F0)^2))
#   }
#   #GB = GenzBretz(maxpts = (p-1)*1e4, abseps = 1e-6, releps = 0)
#   #print(GB$maxpts)
#   #F0 = pmvt(lower = a-mu,upper = b-mu,df = nu,sigma = Sigma)[1]
# 
#   logF0 = pmvt.genz(lower = a-mu,upper = b-mu,nu = nu,sigma = Sigma,uselog2 = TRUE,N = 799)$Estimation
#   F0nnu = pmvt.genz(lower = a-mu,upper = b-mu,nu = nu - 2,sigma = nnu*Sigma,N = 799)$Estimation
# 
#   #Vectors ca and cb
#   SSigma  = nnu*Sigma
#   ssigma2 = diag(SSigma)
#   ca = cb = a0 = a1 = rep(0,p)
#   deltaA  = (nu - 2 + ((a - mu)^2)/diag(SSigma))/(nu - 1)
#   deltaB  = (nu - 2 + ((b - mu)^2)/diag(SSigma))/(nu - 1)
#   yA      = (a - mu)/diag(SSigma)
#   yB      = (b - mu)/diag(SSigma)
# 
#   Wa = Wb = matrix(0,p,p)
# 
#   for(j in 1:p)
#   {
#     #W matrix construction
#       aux1     = onlymeanT0(a = a[-j],b = b[-j],mu = mu[-j] + yA[j]*SSigma[,j][-j],Sigma = deltaA[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),nu = nu-1)
#       ca[j]    = log2prod(dent(a[j],mu[j],ssigma2[j],nu-2),aux1$logF00)
#       Wa[-j,j] = aux1$mean
#       Wa[j,j]  = a[j]
# 
#       aux2     = onlymeanT0(a = a[-j],b = b[-j],mu = mu[-j] + yB[j]*SSigma[,j][-j],Sigma = deltaB[j]*(SSigma[-j,-j] - SSigma[,j][-j]%*%t(SSigma[j,][-j])/ssigma2[j]),nu = nu-1)
#       cb[j]    = log2prod(dent(b[j],mu[j],ssigma2[j],nu-2),aux2$logF00)
#       Wb[-j,j] = aux2$mean
#       Wb[j,j]  = b[j]
#   }
#   muY  = mu + log2ratio(SSigma%*%(ca - cb),logF0)
#   Exx  = muY%*%t(mu) +  log2ratio((F0nnu*diag(p) + Wa%*%diag(ca) - Wb%*%diag(cb))%*%SSigma,logF0)
#   Exx = (Exx + t(Exx))/2
#   varY = Exx - muY%*%t(muY)
#   return(list(mean = muY,EYY = Exx,varcov = varY))
# }
