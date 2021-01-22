########## R function: PoiSimMixModMCMC ##########

# For fitting and inference in a simple (single variance
# parameter) general design Poisson mixed model via 
# mean field variational Bayes.

# Last changed: 31 JUL 2020

PoiSimMixModMCMC <- function(y,X,Z,method,sigmabeta,ssigma,
                             nWarm,nKept,nThin,msgCode)
{
   # Set required dimension variables:

   ncX <- ncol(X)  ;  ncZ <- ncol(Z)  ;  ncC <- ncX + ncZ

   # Form the C matrix:

   Cmat <- cbind(X,Z)

   if (any(method == c("HMC","NUTS")))
   {
      if (msgCode==0) refreshVal <- 0
      if (msgCode==1) refreshVal <- round((nWarm+nKept)/10)
      if (msgCode==2) refreshVal <- round((nWarm+nKept)/100)
      if (msgCode>0)
      {
         # Print message:

         cat("\n")
         if (nThin==1)
         {
            if (method=="HMC")
            {
               cat("   Bayesian density estimation via Hamiltonian Monte Carlo sampling\n")
               cat("   with a warm-up of size",nWarm,"and",nKept,"retained samples.\n")
            }
            if (method=="NUTS")
            {
               cat("   Bayesian density estimation via no-U-turn sampling with\n")
               cat("   a warm-up of size",nWarm,"and",nKept,"retained samples.\n")
            }
         }         
         if (nThin>1)
         {
            if (method=="HMC")
            {
               cat("   Bayesian density estimation via Hamiltonian Monte Carlo sampling\n")
               cat("   with a warm-up of size ",nWarm,", ",nKept," retained samples and a\n",sep="")
               cat("   thinning factor of ",nThin,".\n",sep="")
            }
            if (method=="NUTS")
            {
               cat("   Bayesian density estimation via no-U-turn sampling with\n")
               cat("   a warm-up of size ",nWarm,", ",nKept," retained samples and a\n",sep="")
               cat("   thinning factor of ",nThin,".\n",sep="")
            }
         }   
      }

      ncZ <- ncol(Z)
      allData <- list(n=length(y),ncZ=ncZ,y=y,X=X,Z=Z,sigmabeta=sigmabeta,ssigma=ssigma)
      
      stanObj <- rstan::sampling(stanmodels$PoissonSimpleMixedModel,
                                 data=allData,warmup=nWarm,iter=(nWarm+nKept),
                                 chains=1,thin=nThin,refresh=refreshVal,
                                 algorithm=method,verbose=FALSE)
      
      betauMCMC <- NULL
      for (j in 1:2)
      {
         charVar <- paste("beta[",as.character(j),"]",sep="") 
         betauMCMC <- rbind(betauMCMC,extract(stanObj,charVar,permuted=FALSE))
      }
      for (k in 1:ncZ)
      {
         charVar <- paste("u[",as.character(k),"]",sep="")
         betauMCMC <- rbind(betauMCMC,extract(stanObj,charVar,permuted=FALSE))
      }
      sigmaMCMC <- as.vector(extract(stanObj,"sigma",permuted=FALSE))
   }

   if (method == "slice")
   {
      hyperPars <- c(sigmabeta,ssigma)
      MCMCobj <- PoiSimMixModSlice(y,X,Z,hyperPars,nWarm,nKept,nThin,msgCode)
      betauMCMC <- MCMCobj$betau
      sigmaMCMC <-  MCMCobj$sigma
   }

   return(list(betauMCMC=betauMCMC,sigmaMCMC=sigmaMCMC))
}  

############ End of PoiSimMixModMCMC ############

