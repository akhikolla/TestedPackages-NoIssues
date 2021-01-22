rmvST = function(n,mu = rep(0, length(lambda)),Sigma = diag(length(lambda)),lambda,nu){
  p = length(mu)
  V = rgamma(n,nu/2,nu/2)
  return(matrix(rep(mu,each=n),n,p) + V*rmvSN(n,rep(0,p),Sigma,lambda))
}
