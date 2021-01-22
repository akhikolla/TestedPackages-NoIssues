RcppKmomentN = function(k,a,b,mu,Sigma)
{
  p = length(k)
  
  k = as.matrix(k)
  a = as.matrix(a)
  b = as.matrix(b)
  mu = as.matrix(mu)
  Sigma = as.matrix(Sigma)
  
  out = RcppmomentsN(k,a,b,mu,Sigma)
  mat  = data.frame(cbind(out$index,out$y),row.names = NULL)
  colnames(mat) = c(paste("k",seq(1,p),sep = ""),"E[k]")
  return(mat)
}
