KmomentFESN = function(k,mu,Sigma,lambda,tau)
{
  p = length(mu)
  tautil = tau/sqrt(1+sum(lambda^2))
  if(tautil< -35){
    #print("normal aproximation")
    Delta = sqrtm(Sigma)%*%lambda/sqrt(1+sum(lambda^2))
    return(KmomentFN(k = k,mu = mu - tautil*Delta,Sigma = Sigma - Delta%*%t(Delta)))
  }
  if(all(k == 0)){
    mat  = data.frame(cbind(t(rep(0,p)),1),row.names = NULL)
    colnames(mat) = c(paste("k",seq(1,p),sep = ""),"E[k]")
    return(mat)
  }
  signs = as.matrix(expand.grid(rep(list(c(-1,1)),p)))
  means = matrix(NA,p,2^p)
  ressum  = array(NA,dim = c(prod(k+1),2^p))
  varpsi = lambda/sqrt(1+sum(lambda^2))

  for(i in 1:(2^p)){
    Lambda   = diag(signs[i,])
    mus      = c(Lambda%*%mu)
    Sigmas   = Lambda%*%Sigma%*%Lambda
    SSs      = sqrtm(Sigmas)
    varpsis  = -c(Lambda%*%varpsi)
    Omegas  = cbind(rbind(Sigmas,-t(varpsis)%*%SSs),rbind(-SSs%*%varpsis,1))
    rownames(Omegas) <- colnames(Omegas)

    mom       = as.matrix(Inorm.table(c(k,0),mu = c(mus,tautil),Sigma = Omegas))[,-c(p+1,p+2)]
    ressum[,i] = as.matrix(mom[,-(1:p)])
  }
  mat  = data.frame(cbind(mom[,1:p],matrix(apply(X = ressum,MARGIN = 1,FUN = sum)/pnorm(tautil),prod(k+1),1)),row.names = NULL)
  mat[prod(k+1),p+1] = 1
  colnames(mat) = c(paste("k",seq(1,p),sep = ""),"E[k]")
  return(mat)
}

Inorm.table = function(k,mu,Sigma)
{
  p = length(k)
  if(p==1){
    res  = rev(recintab.RC(k,rep(0,p),rep(Inf,p),mu,Sigma))
    res2 = res
    mat  = data.frame(cbind(seq(k,0),res,res2),row.names = NULL)
    colnames(mat) = c("k","F(k)","E[k]")
    return(mat)
  }
  else{
    t = rep(NA,p)
    for(i in 1:p){
      t[i] = paste("k",i,"=c(",paste(seq(0,k[i]),collapse = ","),")",sep="",collapse = ",")
    }
    res  = recintab.RC(k,rep(0,p),rep(Inf,p),mu,Sigma)
    res2 = res
    eval(parse(text = paste("dimnames(res) = dimnames(res2) = list(",paste(t,collapse = ","),")",sep = "")))

    df1 = as.data.frame.table(res)
    df2 = as.data.frame.table(res2)
    mat1 = data.matrix(df1[prod(k+1):1,], rownames.force = NA)
    mat2 = data.matrix(df2[prod(k+1):1,-(1:p)], rownames.force = NA)
    mat = cbind(mat1,mat2)
    mat[,1:p] = mat[,1:p] - 1
    colnames(mat) = c(paste("k",seq(1,p),sep = ""),"F(k)","E[k]")
    return(mat)
  }
}

####################################################################
recintab.RC = function(kappa,a,b,mu,Sigma)
{
  n = length(kappa)
  seqq = seq_len(n)
  if(n==1)
  {
    M = matrix(data = 0,nrow = kappa+1,ncol = 1)
    s1 = sqrt(Sigma)
    aa = -mu/s1
    M[1] = 1-pnorm(aa)
    if(kappa>0)
    {
      pdfa = s1*dnorm(aa)
      M[2] = mu*M[1]+pdfa
      for(i in seq_len(kappa-1)+1)
      {
        M[i+1] = mu*M[i] + (i-1)*Sigma*M[i-1]
      }
    }
  }else
  {
    #
    #  Create a matrix M, with its nu-th element correpsonds to F_{nu-2}^n(mu,Sigma).
    #
    M  = array(data = 0,dim = kappa+1)
    pk = prod(kappa+1)
    #
    #  We create two long vectors G and H to store the two different sets
    #  of integrals with dimension n-1.
    #
    nn = round(pk/(kappa+1),0)
    begind = cumsum(c(0,nn))
    pk1 = begind[n+1]   # number of (n-1)-dimensional integrals
    #  Each row of cp corresponds to a vector that allows us to map the subscripts
    #  to the index in the long vectors G and H
    cp = matrix(0,n,n)
    for(i in seqq)
    {
      kk = kappa
      kk[i] = 0
      cp[i,] = c(1,cumprod(kk[1:n-1]+1))
    }
    G = rep(0,pk1)
    s = sqrt(diag(Sigma))
    pdfa = dnorm(x = a,mean = mu,sd = s)
    pdfb = dnorm(x = b,mean = mu,sd = s)
    for(i in seqq)
    {
      kappai = kappa[-i]
      ai = a[-i]
      bi = b[-i]
      mui = mu[-i]
      Si = Sigma[-i,i]
      SSi = Sigma[-i,-i] - Si%*%t(Si)/Sigma[i,i]
      ind = (begind[i]+1):begind[i+1]
      if(a[i] != -Inf){
        mai = mui + Si/Sigma[i,i]*(a[i]-mu[i])
        G[ind] = pdfa[i]*recintab.RC(kappai,ai,bi,mai,SSi)
      }
    }
    M[1] = pmvnorm(lower=a-mu,sigma = Sigma)
    cp1 = cp[n,]
    dims  = lapply(X = kappa + 1,FUN = seq_len)
    grid = do.call(expand.grid, dims)

    for(i in seq_len(pk-1)+1)
    {
      kk = as.numeric(grid[i,])
      ii = sum((kk-1)*cp1)+1
      i1 = min(seqq[kk>1])
      kk1 = kk
      kk1[i1] = kk1[i1]-1
      ind3 = ii-cp1[i1]
      M[ii] = mu[i1]*M[ind3]
      for(j in seqq)
      {
        kk2 = kk1[j]-1
        if(kk2 > 0)
        {
          M[ii] = M[ii]+Sigma[i1,j]*kk2*M[ind3-cp1[j]]
        }
        ind4 = as.numeric(begind[j] + sum(cp[j,]*(kk1-1)) - cp[j,j]*kk2+1)
        M[ii] = M[ii]+Sigma[i1,j]*a[j]^kk2*G[ind4]
      }
    }
  }
  return(M)
}

# p = 5
# a  = seq(-0.9,0,length.out = p)
# b  = seq(0,0.9,length.out = p) + seq(0.1,0.6,length.out = p)
# mu = seq(0.25,0.75,length.out = p)*b
# s  = matrix(1.2*rnorm(p^2),p,p)
# Sigma = S = s%*%t(s)
# lambda = 1:p - round((p+1)/2)
# tau=1
#
#
# microbenchmark(
#   momentsFESN(kappa = c(2,1,0,0,1),mu,Sigma,lambda,tau),
#   KmomentFESN(k = c(2,1,0,0,1),mu,Sigma,lambda,tau),times = 2)


# momentsFESN = function(kappa,mu,Sigma,lambda,tau){
#   p = length(mu)
#   if(p==1){
#     means = recintabESN(kappa,rep(0,p),rep(Inf,p),mu,Sigma,lambda,tau) + recintabESN(kappa,rep(0,p),rep(Inf,p),-mu,Sigma,-lambda,tau)
#   }else{
#     signs = as.matrix(expand.grid(rep(list(c(-1,1)),p)))
#     means = matrix(0,prod(kappa+1),p+1)
#     for(i in 1:(2^p))
#     {
#       Lambda   = diag(signs[i,])
#       mus      = c(Lambda%*%mu)
#       Ss       = Lambda%*%Sigma%*%Lambda
#       lambdas  = c(Lambda%*%lambda)
#       res      = recintabESN(kappa,rep(0,p),rep(Inf,p),mus,Ss,lambdas,tau)
#       means  = means + res
#     }
#   }
#   means[1,p+1] = 1
#   Moments = means[,p+1]
#   return(cbind(means[,1:p]/2^p,Moments))
# }
#
# recintabESN0 = function(kappa,a,b,mu,Sigma,lambda,tau)
# {
#   p = length(kappa)
#   seqq = seq_len(p)
#   if(p==1)
#   {
#     M = matrix(data = NA,nrow = kappa+1,ncol = 1)
#     M[1] = AcumESN(b,mu,Sigma,lambda,tau) - AcumESN(a,mu,Sigma,lambda,tau)
#     if(kappa>0)
#     {
#       #normal moments
#       s = sqrt(Sigma)
#       tautil = tau/sqrt(1+sum(lambda^2))
#       phi    = lambda/sqrt(1+sum(lambda^2))
#       eta    = dnorm(tau)/pnorm(tautil)/sqrt(1+sum(lambda^2))*exp(phi^2*tau^2/2)
#       Gamma  = Sigma/(1+lambda^2)
#       mub    = lambda*tau*Gamma/s
#       N = recintab0(kappa-1,a,b,mu-mub,Gamma)
#
#       #ESN moments
#       da     = dmvESN1(a,mu,Sigma,lambda,tau)
#       db     = dmvESN1(b,mu,Sigma,lambda,tau)
#       normct = lambda*s*eta
#       M[2] = mu*M[1] + Sigma*(da-db) + normct*N[1]
#       if(a == -Inf){a = 0}
#       if(b == Inf){b = 0}
#       for(i in seq_len(kappa-1)+1)
#       {
#         da = da*a
#         db = db*b
#         M[i+1] = mu*M[i] + (i-1)*Sigma*M[i-1] + Sigma*(da - db) + normct*N[i]
#       }
#     }
#     return(M)
#   }else
#   {
#     #
#     #  Create a matrix M, with its nu-th element correpsonds to F_{nu-2}^n(mu,Sigma).
#     #
#     #M  = array(data = 0,dim = kappa+1)
#     pk = prod(kappa+1)
#     M = rep(0,pk)
#     M[1] = pmvESN(a,b,mu,Sigma,lambda,tau)
#     if(sum(kappa) == 0){return(M)}
#     #
#     #  We create two long vectors G and H to store the two different sets
#     #  of integrals with dimension p-1.
#     #
#     nn = round(pk/(kappa+1),0)
#     begind = cumsum(c(0,nn))
#     pk1 = begind[p+1]   # number of (p-1)-dimensional integrals
#     #  Each row of cp corresponds to a vector that allows us to map the subscripts
#     #  to the index in the long vectors G and H
#     cp = matrix(0,p,p)
#     pdfa = pdfb = rep(NA,p)
#
#     for(i in seqq)
#     {
#       kk = kappa
#       kk[i] = 0
#       cp[i,] = c(1,cumprod(kk[1:p-1]+1))
#     }
#     G = rep(0,pk1)
#     H = rep(0,pk1)
#     s = sqrt(diag(Sigma))                            # usamos s?
#     a1   = a-mu
#     b1   = b-mu
#     SS   = sqrtm(Sigma)
#     iSS  = solve(SS)
#
#     #aux
#     tautil = tau/sqrt(1+sum(lambda^2))
#     Phi    = diag(p) + lambda%*%t(lambda)
#     Gamma  = SS%*%solve(Phi)%*%SS
#     iGamma = iSS%*%Phi%*%iSS
#     varphi = iSS%*%lambda
#     mub    = tau*Gamma%*%varphi
#     eta    = dnorm(tau)/pnorm(tautil)/sqrt(1+sum(lambda^2))*as.numeric(exp(t(mub)%*%iGamma%*%mub/2))
#     delta  = eta*Sigma%*%varphi
#
#     for(i in seqq)
#     {
#       SSi = Sigma[-i,-i] - Sigma[-i,i]%*%t(Sigma[-i,i])/Sigma[i,i]
#       cj = (1 + t(varphi[-i])%*%SSi%*%varphi[-i])^{-1/2}
#       tilvarphi = varphi[i] + (Sigma[-i,i]/Sigma[i,i])%*%varphi[-i]
#       #pdfs in d_kappa vector
#       pdfa[i] = dmvESN1(y = a1[i],Sigma = Sigma[i,i],lambda = sqrt(Sigma[i,i])*cj*tilvarphi,tau = cj*tau)
#       pdfb[i] = dmvESN1(y = b1[i],Sigma = Sigma[i,i],lambda = sqrt(Sigma[i,i])*cj*tilvarphi,tau = cj*tau)
#       #more
#       ind = (begind[i]+1):begind[i+1]
#       if(a[i] != -Inf){
#         mai = mu[-i] + Sigma[-i,i]/Sigma[i,i]*a1[i]
#         G[ind] = pdfa[i]*recintabESN0(kappa[-i],a[-i],b[-i],mai,SSi,
#                                       c(sqrtm(SSi)%*%varphi[-i]),c(tau + tilvarphi*a1[i]))
#       }
#       if(b[i] != Inf){
#         mbi = mu[-i] + Sigma[-i,i]/Sigma[i,i]*b1[i]
#         H[ind] = pdfb[i]*recintabESN0(kappa[-i],a[-i],b[-i],mbi,SSi,
#                                       c(sqrtm(SSi)%*%varphi[-i]),c(tau + tilvarphi*b1[i]))
#       }
#     }
#
#     N = recintab0(kappa,a,b,c(mu-mub),Gamma)
#
#     a[a == -Inf] = 0
#     b[b == Inf] = 0
#     cp1 = cp[p,]
#     dims  = lapply(X = kappa + 1,FUN = seq_len)
#     grid = do.call(expand.grid, dims)
#
#     for(i in seq_len(pk-1)+1)
#     {
#       kk = as.numeric(grid[i,])
#       ii = sum((kk-1)*cp1)+1
#       i1 = min(seqq[kk>1])
#       kk1 = kk
#       kk1[i1] = kk1[i1]-1
#       ind3 = ii-cp1[i1]
#       M[ii] = mu[i1]*M[ind3] + delta[i1]*N[ind3]  #here we add the normal part
#       for(j in seqq)
#       {
#         kk2 = kk1[j]-1
#         if(kk2 > 0)
#         {
#           M[ii] = M[ii]+Sigma[i1,j]*kk2*M[ind3-cp1[j]]
#         }
#         ind4 = as.numeric(begind[j] + sum(cp[j,]*(kk1-1)) - cp[j,j]*kk2+1)
#         M[ii] = M[ii]+Sigma[i1,j]*(a[j]^kk2*G[ind4]-b[j]^kk2*H[ind4])
#       }
#     }
#     return(M)
#   }
# }
#
# #############################################################################
#
# recintabESN = function(kappa,a,b,mu,Sigma,lambda,tau)
# {
#   p = length(kappa)
#   seqq = seq_len(p)
#   if(p==1)
#   {
#     M = matrix(data = NA,nrow = kappa+1,ncol = 1)
#     M[1] = AcumESN(b,mu,Sigma,lambda,tau) - AcumESN(a,mu,Sigma,lambda,tau)
#     if(kappa>0)
#     {
#       #normal moments
#       s = sqrt(Sigma)
#       tautil = tau/sqrt(1+sum(lambda^2))
#       phi    = lambda/sqrt(1+sum(lambda^2))
#       eta    = dnorm(tau)/pnorm(tautil)/sqrt(1+sum(lambda^2))*exp(phi^2*tau^2/2)
#       Gamma  = Sigma/(1+lambda^2)
#       mub    = lambda*tau*Gamma/s
#       N = recintab0(kappa-1,a,b,mu-mub,Gamma)
#
#       #ESN moments
#       da     = dmvESN1(a,mu,Sigma,lambda,tau)
#       db     = dmvESN1(b,mu,Sigma,lambda,tau)
#       normct = lambda*s*eta
#       M[2] = mu*M[1] + Sigma*(da-db) + normct*N[1]
#       if(a == -Inf){a = 0}
#       if(b == Inf){b = 0}
#       for(i in seq_len(kappa-1)+1)
#       {
#         da = da*a
#         db = db*b
#         M[i+1] = mu*M[i] + (i-1)*Sigma*M[i-1] + Sigma*(da - db) + normct*N[i]
#       }
#     }
#     return(cbind(seq_len(kappa+1)-1,M))
#   }else
#   {
#     #
#     #  Create a matrix M, with its nu-th element correpsonds to F_{nu-2}^n(mu,Sigma).
#     #
#     #M  = array(data = 0,dim = kappa+1)
#     pk = prod(kappa+1)
#     M = rep(0,pk)
#     M[1] = pmvESN(a,b,mu,Sigma,lambda,tau)
#     if(sum(kappa) == 0){return(cbind(matrix(0,1,p),M))}
#     #
#     #  We create two long vectors G and H to store the two different sets
#     #  of integrals with dimension p-1.
#     #
#     nn = round(pk/(kappa+1),0)
#     begind = cumsum(c(0,nn))
#     pk1 = begind[p+1]   # number of (p-1)-dimensional integrals
#     #  Each row of cp corresponds to a vector that allows us to map the subscripts
#     #  to the index in the long vectors G and H
#     cp = matrix(0,p,p)
#     pdfa = pdfb = rep(NA,p)
#
#     for(i in seqq)
#     {
#       kk = kappa
#       kk[i] = 0
#       cp[i,] = c(1,cumprod(kk[1:p-1]+1))
#     }
#     G = rep(0,pk1)
#     H = rep(0,pk1)
#     s = sqrt(diag(Sigma))                            # usamos s?
#     a1   = a-mu
#     b1   = b-mu
#     SS   = sqrtm(Sigma)
#     iSS  = solve(SS)
#
#     #aux
#     tautil = tau/sqrt(1+sum(lambda^2))
#     Phi    = diag(p) + lambda%*%t(lambda)
#     Gamma  = SS%*%solve(Phi)%*%SS
#     iGamma = iSS%*%Phi%*%iSS
#     varphi = iSS%*%lambda
#     mub    = tau*Gamma%*%varphi
#     eta    = dnorm(tau)/pnorm(tautil)/sqrt(1+sum(lambda^2))*as.numeric(exp(t(mub)%*%iGamma%*%mub/2))
#     delta  = eta*Sigma%*%varphi
#
#     for(i in seqq)
#     {
#       SSi = Sigma[-i,-i] - Sigma[-i,i]%*%t(Sigma[-i,i])/Sigma[i,i]
#       cj = (1 + t(varphi[-i])%*%SSi%*%varphi[-i])^{-1/2}
#       tilvarphi = varphi[i] + (Sigma[-i,i]/Sigma[i,i])%*%varphi[-i]
#       #pdfs in d_kappa vector
#       pdfa[i] = dmvESN1(y = a1[i],Sigma = Sigma[i,i],lambda = sqrt(Sigma[i,i])*cj*tilvarphi,tau = cj*tau)
#       pdfb[i] = dmvESN1(y = b1[i],Sigma = Sigma[i,i],lambda = sqrt(Sigma[i,i])*cj*tilvarphi,tau = cj*tau)
#       #more
#       ind = (begind[i]+1):begind[i+1]
#       if(a[i] != -Inf){
#         mai = mu[-i] + Sigma[-i,i]/Sigma[i,i]*a1[i]
#         G[ind] = pdfa[i]*recintabESN0(kappa[-i],a[-i],b[-i],mai,SSi,
#                                       c(sqrtm(SSi)%*%varphi[-i]),c(tau + tilvarphi*a1[i]))
#       }
#       if(b[i] != Inf){
#         mbi = mu[-i] + Sigma[-i,i]/Sigma[i,i]*b1[i]
#         H[ind] = pdfb[i]*recintabESN0(kappa[-i],a[-i],b[-i],mbi,SSi,
#                                       c(sqrtm(SSi)%*%varphi[-i]),c(tau + tilvarphi*b1[i]))
#       }
#     }
#
#     N = recintab0(kappa,a,b,c(mu-mub),Gamma)
#
#     a[a == -Inf] = 0
#     b[b == Inf] = 0
#     cp1 = cp[p,]
#     dims  = lapply(X = kappa + 1,FUN = seq_len)
#     grid = do.call(expand.grid, dims)
#
#     for(i in seq_len(pk-1)+1)
#     {
#       kk = as.numeric(grid[i,])
#       ii = sum((kk-1)*cp1)+1
#       i1 = min(seqq[kk>1])
#       kk1 = kk
#       kk1[i1] = kk1[i1]-1
#       ind3 = ii-cp1[i1]
#       M[ii] = mu[i1]*M[ind3] + delta[i1]*N[ind3]  #here we add the normal part
#       for(j in seqq)
#       {
#         kk2 = kk1[j]-1
#         if(kk2 > 0)
#         {
#           M[ii] = M[ii]+Sigma[i1,j]*kk2*M[ind3-cp1[j]]
#         }
#         ind4 = as.numeric(begind[j] + sum(cp[j,]*(kk1-1)) - cp[j,j]*kk2+1)
#         M[ii] = M[ii]+Sigma[i1,j]*(a[j]^kk2*G[ind4]-b[j]^kk2*H[ind4])
#       }
#     }
#   }
#   return(cbind(grid-1,round(M,4)))
# }
#
#
# #####################################################
# #####################################################
