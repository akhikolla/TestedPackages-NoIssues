########## R function: checkChains ##########

# For checking the chains of a stochastic densEstBayes() 
# fit object:

# Last changed: 14 SEP 2020

checkChains <- function(fitObject,locations="equally-spaced")
{
   # Check legality of the "locations" argument:

   if (!any(locations==c("equally-spaced","hexiles")))
   {
      warnStr1 <- "The inputted locations value  was neither \"equally-spaced\" nor \"hexiles\"."
      warnStr2 <- "The default value of \"equally-spaced\" was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      locations <- "equally-spaced"
   }

   # Extract type: 

   method <- fitObject$method
   if (method=="SMFVB")
      stop("There are no chains to check when method=\"SMFVB\".\n")
   
   # Extract fit information:

   stochaFitObj <- fitObject$stochaFitObj
   range.x <- fitObject$range.x
   intKnots <- fitObject$intKnots
   betauMCMC <- stochaFitObj$betauMCMC

   # Set grid-wise matrices: 

   gridSize <- 1001
   anchorLow <- -0.05  ;  anchorUpp <- 1.05
   xTrang <- seq(anchorLow,anchorUpp,length=gridSize)
   Xg <- cbind(1,xTrang)
   Zg <- ZOSull(xTrang,intKnots=intKnots,range.x=c(anchorLow,anchorUpp))
   Cg <- cbind(Xg,Zg)
   etagMCMC <- crossprod(t(Cg),betauMCMC)

   if (locations=="equally-spaced")
   {
      # Form chains corresponding to equally-spaced vertical slices:

      equalSpaceGrid <- seq(anchorLow,anchorUpp,length=7)[-c(1,7)]
      igVec <- rep(NA,5)
      for (iq in 1:5)
         igVec[iq] <- which.min((xTrang-equalSpaceGrid[iq])^2)
      MCMCmat <- etagMCMC[igVec[1],]
      for (j in 2:5)
         MCMCmat <- cbind(MCMCmat,etagMCMC[igVec[j],])
   }

   if (locations=="hexiles")
   {
      # Form chains corresponding to vertical slices at the sample hexiles:

      igVec <- rep(NA,5)
      for (iq in 1:5)
         igVec[iq] <- which.min((xTrang-fitObject$sampHexTran[iq])^2)
      MCMCmat <- etagMCMC[igVec[1],]
      for (j in 2:5)
         MCMCmat <- cbind(MCMCmat,etagMCMC[igVec[j],])
   }

   # Produce graphical summary:

   if (locations=="equally-spaced")
      parNamesVal <- c("1/6 of range","2/6 of range","3/6 of range",
                       "4/6 of range","5/6 of range")
   if (locations=="hexiles")
      parNamesVal <- c("1st hexile","2nd hexile","3rd hexile",
                       "4th hexile","5th hexile")
   summChainsDensEst(MCMCmat,parNames=parNamesVal)

   invisible()
}

############ End of checkChains ############
