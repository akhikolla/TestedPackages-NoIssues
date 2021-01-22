smoothrhos <- function(rho, BasiseBsc, zero_range){ 
    
    neBsc = length(rho)
    x     = seq(1,neBsc,length = neBsc)
    PP    = BasiseBsc$eigenvectorsQR
    ee    = as.vector(neBsc * BasiseBsc$eigenvalues)
    coefs = t(PP) %*% rho

    temp  = get_lambda(zero_range = zero_range, coefs = coefs, eigens = ee, rhos = 1 , neBsc = neBsc, q = 2)         
    rho.fit = PP %*% (coefs / (1+(10^temp) * ee)) 
    pmax( 10^-5, rho.fit)    
  }
