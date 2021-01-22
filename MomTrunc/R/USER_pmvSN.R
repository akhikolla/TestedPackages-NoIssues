#######################################################################################
#######################################################################################

pmvSN = function(lower = rep(-Inf,length(lambda)),upper=rep(Inf,length(lambda)),mu = rep(0,length(lambda)),Sigma,lambda,log2 = FALSE){
  return(pmvESN(lower = lower,upper = upper,mu = mu,Sigma = Sigma,lambda = lambda,tau = 0,log2 = log2))
}
