#######################################################################################
#######################################################################################
dmvST = function(x,mu=rep(0,length(lambda)),Sigma=diag(length(lambda)),lambda,nu){
  return(dmvEST(x,mu,Sigma,lambda,tau=0,nu))
}
