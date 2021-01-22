#Bivariate T density

dt2d = function(pair, rho = 0, nu = 4){
  if(any(is.infinite(pair))){return(0)}
  x = pair[1]
  y = pair[2]
  if (nu == Inf) return(dmvnorm(x = c(x,y),sigma = matrix(c(1,rho,rho,1),2,2))) 
  
  # Argument:
  xoy <- (x^2 - 2*rho*x*y + y^2)/ (2*(1 - rho^2))
  
  # Density:
  density <- (1 + 2*xoy/nu)^(-(nu+2)/2) / (2*pi*sqrt(1-rho^2))
  # Return value:
  density
}
