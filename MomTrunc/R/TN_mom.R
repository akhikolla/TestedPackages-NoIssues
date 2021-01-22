KmomentN = function(k,a,b,mu,Sigma)
{
  p = length(k)
  if(p==1){
    res  = rev(recintab(k,a,b,mu,Sigma))
    res2 = res/as.numeric(res)[1]
    mat  = data.frame(cbind(seq(k,0),res,res2),row.names = NULL)
    colnames(mat) = c("k","F(k)","E[k]")
    return(mat)
  }
  else{
    t = rep(NA,p)
    for(i in 1:p){
      t[i] = paste("k",i,"=c(",paste(seq(0,k[i]),collapse = ","),")",sep="",collapse = ",")
    }
    res  = recintab(k,a,b,mu,Sigma)
    res2 = res/as.numeric(res)[1]
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

######################################################

recintab = function(kappa,a,b,mu,Sigma)
{
  n = length(kappa)
  seqq = seq_len(n)
  if(n==1)
  {
    M = matrix(data = 0,nrow = kappa+1,ncol = 1)
    s1 = sqrt(Sigma)
    aa = (a-mu)/s1
    bb = (b-mu)/s1
    M[1] = pnorm2(bb)-pnorm(aa)
    if(kappa>0)
    {
      pdfa = s1*dnorm(aa)
      pdfb = s1*dnorm(bb)
      M[2] = mu*M[1]+pdfa-pdfb
      if(a == -Inf){a = 0}
      if(b == Inf){b = 0}
      for(i in seq_len(kappa-1)+1)
      {
        pdfa = pdfa*a
        pdfb = pdfb*b
        M[i+1] = mu*M[i] + (i-1)*Sigma*M[i-1] + pdfa - pdfb
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
    H = rep(0,pk1)
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
        G[ind] = pdfa[i]*recintab(kappai,ai,bi,mai,SSi)
      }
      if(b[i] != Inf){
        mbi = mui+Si/Sigma[i,i]*(b[i]-mu[i])
        H[ind] = pdfb[i]*recintab(kappai,ai,bi,mbi,SSi)
      }
    }
    M[1] = pmvnormt(lower=a, upper=b, mean=mu,sigma = Sigma)
    
    a[a == -Inf] = 0
    b[b == Inf] = 0
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
        M[ii] = M[ii]+Sigma[i1,j]*(a[j]^kk2*G[ind4]-b[j]^kk2*H[ind4])
      }
    }
  }
  return(M)
}


