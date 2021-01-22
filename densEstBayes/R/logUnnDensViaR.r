########## R function: logUnnDensViaR ##########

# Evaluation of the logarithm of the unnormalised
# density function of the current `betau' parameter.

# Last changed: 21 JUL 2020

logUnnDensViaR <- function(j,betau.j,betau.not.j,Cmat,CTy,sigsq.betau)
{
   bFun <- function(x) return(exp(x))
   ncC <- ncol(Cmat)                            
   lenbetau <- 1 + length(betau.not.j)

   # Do first term:

   firstTerm <- CTy[j]*betau.j

   # Do second term:

   betauOutOrder <- c(betau.j,betau.not.j)       
   if (j==1) idx <- 1:ncC
   if ((j>1)&(j<ncC)) idx <- c((2:j),1,(j+1):ncC)
   if (j==ncC) idx <- c(2:ncC,1)
   betau <- betauOutOrder[idx]
   
   secondTerm <- -sum(bFun(as.vector(crossprod(t(Cmat),betau))))

   # Do third term:
   
   thirdTerm <- -betau.j^2/(2*sigsq.betau[j])

   # Return answer:

   return(firstTerm + secondTerm + thirdTerm)
}     

############ End of logUnnDensViaR ############

