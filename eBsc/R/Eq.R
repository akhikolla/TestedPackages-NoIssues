Eq<-function(log10lambda,coefs,eigens,rhos,neBsc,q) 
  {
    lambda = 10^log10lambda
    r.mat1 = 2 * q * log((eigens)^( 1 / (2 * q)))/(1 + eigens * lambda * rhos)
    r.mat = eigens * lambda / (1 + eigens * lambda * rhos)
    sum((coefs ^ 2 * r.mat * r.mat1)[-c(1:q)]) - sum(r.mat1[-c(1:q)]) * (sum(coefs^2 * r.mat)+1) / (neBsc+1)
}

