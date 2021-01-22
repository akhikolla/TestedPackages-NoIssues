RcppKmomentT = function(k,a,b,mu,Sigma,nu)
{
  p = length(mu)
  if(k == 0){
    mat  = data.frame(matrix(c(rep(0,p),1),1,p+1),row.names = NULL)
    colnames(mat) = c(paste("k",seq(1,p),sep = ""),"E[k]")
    return(mat)
  }else{
    out = momentsT(k,a,b,mu,Sigma,nu)
    mat = Rcpporder(k,cbind(out$index,out$y))
    mat  = data.frame(mat,row.names = NULL)
    colnames(mat) = c(paste("k",seq(1,p),sep = ""),"E[k]")
    return(mat)
    }
}

#######################################################%

momentsT = function(kappa,a,b,mu,S,nu){
  out = hfunrec2(kappa,a,b,mu,nu*S,nu) 
  out$y = out$y/out$y[1]
  return(out)
}
