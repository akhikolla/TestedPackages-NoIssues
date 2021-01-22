KmomentT = function(k,a,b,mu,Sigma,nu)
{
  p = length(k)
  if(all(k == 0))
  {
    res  = pmvt.genz(lower = a - mu,upper = b - mu,nu = nu,sigma = Sigma)$Estimation
    mat  = data.frame(cbind(rep(0,p),round(res,4),1),row.names = NULL)
    colnames(mat) = c(paste("k",seq(1,p),sep = ""),"F(k)","E[k]")
    return(mat)
  }else{
    res  = Fk(k,a,b,mu,Sigma,nu)
    mat  = cbind(res$ks,c(res$ress))
    F0   = mat[dim(mat)[1],dim(mat)[2]]
    mat  = data.frame(cbind(round(mat,4),round(mat[,dim(mat)[2]]/F0,4)),row.names = NULL)
    colnames(mat) = c(paste("k",seq(1,p),sep = ""),"F(k)","E[k]")
    return(mat)}
}
