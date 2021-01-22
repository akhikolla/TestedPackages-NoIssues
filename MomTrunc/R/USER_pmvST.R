#######################################################################################
#######################################################################################

pmvST = function(lower = rep(-Inf,length(lambda)),upper=rep(Inf,length(lambda)),mu = rep(0,length(lambda)),Sigma,lambda,nu,log2 = FALSE){
  return(pmvEST(lower = lower,upper = upper,mu = mu,Sigma = Sigma,lambda = lambda,tau = 0,nu = nu,log2 = log2))
}
