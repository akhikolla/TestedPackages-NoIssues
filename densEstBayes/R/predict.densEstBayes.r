########## R function: predict.densEstBayes ##########

# For prediction for new data from a densEstBayes() fit object:

# Last changed: 14 SEP 2020

predict.densEstBayes <- function(object,newdata,cred.fit=FALSE,credLev=0.95,...)
{
   fitObject <- object
   xNew <- newdata

   # Check legality of the credLev parameter:

   if (credLev<=0) 
   {
      warnStr1 <- "The inputted credible level is zero or negative."
      warnStr2 <- "The default value of 0.95 was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      credLev <- 0.95
   }
   if (credLev>=1) 
   {
      warnStr1 <- "The inputted credible level is 1 or higher."
      warnStr2 <- "The default value of 0.95 was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      credLev <- 0.95
   }

   # Extract type: 

   method <- fitObject$method

   # Extract fit values:

   determFitObj <- fitObject$determFitObj
   stochaFitObj <- fitObject$stochaFitObj

   # Extract spline basis function parameters:

   range.x <- fitObject$range.x
   intKnots <- fitObject$intKnots

   if (method=="SMFVB")
   {
      mu.q.betau <- determFitObj$mu.q.betau
      Sigma.q.betau <- determFitObj$Sigma.q.betau
   }
   if (any(method==c("HMC","NUTS","slice")))
       betauMCMC <- stochaFitObj$betauMCMC

   # First the estimate over the plot.densEstBayes() plotting
   # grid is computed due to the need for normalisation:

   anchorLow <- -0.05  ;  anchorUpp <- 1.05
   xTrang <- seq(anchorLow,anchorUpp,length=1001)
   Xg <- cbind(1,xTrang)
   Zg <- ZOSull(xTrang,intKnots=intKnots,range.x=c(anchorLow,anchorUpp))
   Cg <- cbind(Xg,Zg)

   if (method=="SMFVB")
   {
      etaHatg <- as.vector(crossprod(t(Cg),mu.q.betau))
      sdg <- sqrt(rowSums(crossprod(t(Cg),Sigma.q.betau)*Cg))
      lowg <- etaHatg - qnorm((1+credLev)/2)*sdg
      uppg <- etaHatg + qnorm((1+credLev)/2)*sdg
   } 

   if (any(method==c("HMC","NUTS","slice")))
   {
      etaHatMCMC <- crossprod(t(Cg),betauMCMC)
      etaHatg <- apply(etaHatMCMC,1,mean)
      lowg <- apply(etaHatMCMC,1,quantile,(1-credLev)/2)
      uppg <- apply(etaHatMCMC,1,quantile,(1+credLev)/2)
   } 

   densEstUnng <- exp(etaHatg)
   normFac <- trapint(xTrang,densEstUnng)

   densEstTrang <- densEstUnng/normFac
   densLowTrang <- exp(lowg)/normFac
   densUppTrang <- exp(uppg)/normFac
  
   # Transform `xg' and 'densEstg' to the original units:

   xg <- range.x[1] + (range.x[2] - range.x[1])*xTrang
   densEstg <- densEstTrang/(range.x[2] - range.x[1])
   densLowg <- densLowTrang/(range.x[2] - range.x[1])
   densUppg <- densUppTrang/(range.x[2] - range.x[1])
  
   # Now use linear interpolation to obtain predictions at
   # the new data:

   # Check the legality of the "newdata" input:

   if (any(xNew<min(xg))) stop("At least one abscissae value is below default plotting grid.\n")
   if (any(xNew>max(xg))) stop("At least one abscissae value is above default plotting grid.\n")

   predictor <- approx(xg,densEstg,xNew)$y
   credLow <- approx(xg,densLowg,xNew)$y
   credUpp <- approx(xg,densUppg,xNew)$y

   return(list(fit = predictor,credLow.fit=credLow,credUpp.fit=credUpp))
}

############ End of predict.densEstBayes ############
