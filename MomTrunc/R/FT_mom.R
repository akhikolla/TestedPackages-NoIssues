KmomentFT = function(k,mu,Sigma,nu)
{
  p = length(mu)
  if(all(k == 0)){
    mat  = data.frame(cbind(t(rep(0,p)),1),row.names = NULL)
    colnames(mat) = c(paste("k",seq(1,p),sep = ""),"E[k]")
    return(mat)
  }
  if(p==1)
  {
    mat  = data.frame(cbind(seq(k,0),as.matrix(KmomentT(k,a = 0,b = 10^10,-mu,Sigma,nu))[,2]+as.matrix(KmomentT(k,a = 0,b = 10^10,mu,Sigma,nu))[,2]),row.names = NULL)
    colnames(mat) = c("k","E[k]")
    return(mat)

  }else{
    signs = as.matrix(expand.grid(rep(list(c(-1,1)),p)))
    means = matrix(NA,p,2^p)
    ressum  = array(NA,dim = c(sum(k)+1,2^p))
    for(i in 1:(2^p))
    {
      Lambda   = diag(signs[i,])
      mus      = c(Lambda%*%mu)
      Ss       = Lambda%*%Sigma%*%Lambda
      mom       = as.matrix(KmomentT(k,a = rep(0,p),b = rep(10^10,p),mus,Ss,nu))[,-(p+2)]
      ressum[,i] = as.matrix(mom[,-(1:p)])
    }
    mat  = data.frame(cbind(mom[,1:p],matrix(apply(X = ressum,MARGIN = 1,FUN = sum),sum(k)+1,1)),row.names = NULL)
    mat[sum(k)+1,p+1] = 1
    colnames(mat) = c(paste("k",seq(1,p),sep = ""),"E[k]")
    return(mat)
  }
}
