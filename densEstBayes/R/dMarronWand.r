########## R-function: dMarronWand ##########

# For computing the ordinate values of a
# Marron & Wand normal mixture density,
# distribution function or derivatives.

# Last changed: 23 OCT 2018

dMarronWand <- function(x,densNum,drv=0)
{
   parms <- MarronWandParm(densNum)
   w <- parms$w
   mu <- parms$mu
   sigmasq <- parms$sigmasq
  
   k <- length(w)
   sig <- sqrt(sigmasq)
   ans <- 0
   for (j in 1:k) 
   {
      if (drv==0)    
         newTerm <- w[j]*dnorm(x,mu[j],sig[j])
      if (drv==(-1))    
         newTerm <- w[j]*pnorm(x,mu[j],sig[j])
      if (drv==1) 
         newTerm <- (-w[j]*((x-mu[j])/sig[j])
                     *dnorm(x,mu[j],sig[j])/sig[j])
      if (drv==2) 
         newTerm <- (w[j]*(((x-mu[j])/sig[j])^2-1)
                     *dnorm(x,mu[j],sig[j])/(sig[j]^2))
      
      ans <- ans + newTerm
   }
   return(ans)
} 

######## End S-function dMarronWand ##########


