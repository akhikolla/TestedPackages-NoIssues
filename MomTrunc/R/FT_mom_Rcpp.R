RcppKmomentFT = function(k,mu,Sigma,nu)
{
  p = length(mu)
  if(k == 0){
    mat  = data.frame(matrix(c(rep(0,p),1),1,p+1),row.names = NULL)
    colnames(mat) = c(paste("k",seq(1,p),sep = ""),"E[k]")
    return(mat)
  }else{
    out = RcppmomentsFT(k,mu,Sigma,nu)
    out  = data.frame(out,row.names = NULL)
    colnames(out) = c(paste("k",seq(1,p),sep = ""),"E[k]")
    return(out)
  }
}
