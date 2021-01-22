RcppKmomentFN = function(k,mu,Sigma)
{
  p = length(mu)
  if(all(k == 0)){
    mat  = data.frame(matrix(c(rep(0,p),1),1,p+1),row.names = NULL)
    colnames(mat) = c(paste("k",seq(1,p),sep = ""),"E[k]")
    return(mat)
  }else{
    out = RcppmomentsFN(k,mu,Sigma)
    out  = data.frame(cbind(out$index,out$y),row.names = NULL)
    colnames(out) = c(paste("k",seq(1,p),sep = ""),"E[k]")
    return(out)
  }
}
