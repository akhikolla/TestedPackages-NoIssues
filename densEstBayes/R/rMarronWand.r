########## R-function: rMarronWand ##########

# For computing a random sample from a 
# Marron & Wand Normal Mixture density.

# Last changed: 23 OCT 2018

rMarronWand <- function(n,densNum)
{
   # Check legality of inputs:

   n <- round(n)
   if (n<=0) stop("n must be a positive integer.\n")
   
   densNum <- round(densNum)
   if (!any(densNum==(1:15))) 
      stop("densNum must be an integer between 1 and 15 inclusive.\n")

   parms <- MarronWandParm(densNum)
   w <- parms$w
   mu <- parms$mu
   sigmasq <- parms$sigmasq
  
   return(nmsamp(n,w,mu,sigmasq))
} 

######## End S-function rMarronWand ##########


