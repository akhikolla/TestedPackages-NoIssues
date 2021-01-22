########## R function: PoiSimMixModSlice ##########

# For performing the Markov chain Monte Carlo
# iterations, based on slice sampling, for 
# a simple Poisson mixed model.

# Last changed: 27 JUL 2020

PoiSimMixModSlice <- function(y,X,Z,hyperPars,nWarm,nKept,nThin,msgCode)
{
   # Set MCMC dimension variable: 
  
   numMCMC <-  nWarm + nKept 

   # Set hyperpameters:

   sigmabetaHYP <- hyperPars[1]
   AHYP <- hyperPars[2]

   # Set dimension and constant matrices:

   ncX <- ncol(X)
   Cmat <- cbind(X,Z)
   CTy <- crossprod(Cmat,y)

   # If 'msgCode' is positive then print an informational message:

   if (msgCode>0)
   {
      cat("\n")
      if (nThin==1)
      {
         cat("   Bayesian density estimation via slice sampling with\n")
         cat("   a warm-up of size ",nWarm," and ",nKept," retained samples.\n",sep="")
      }
      if (nThin>1)
      {
         cat("   Bayesian density estimation via slice sampling with\n")
         cat("   a warm-up of size ",nWarm,", ",nKept," retained samples and a\n",sep="")
         cat("   thinning factor of ",nThin,".\n",sep="")
      }
   }

   # Obtain slice-based MCMC samples:

   innerObj <- PoiSMMsliceInner(numMCMC,ncX,y,Cmat,CTy,
                                sigmabetaHYP,AHYP,msgCode)
                                      
   # Extract samples and discard the warm-up values:

   betauMCMC <- innerObj$betau[,-(1:(nWarm))]
   sigmaMCMC <- sqrt(innerObj$sigsq)[-(1:(nWarm))]
 
   # Do thinning if 'nThin' exceeds unity:

   if (nThin>1)
   {
      thinnedInds <- seq(1,nKept,by=nThin)
      betauMCMC <- betauMCMC[,thinnedInds]
      sigmaMCMC <- sigmaMCMC[thinnedInds]
   }

   # Return kept samples:

   return(list(betau=betauMCMC,sigma=sigmaMCMC))
}

############ End of PoiSimMixModSlice ############