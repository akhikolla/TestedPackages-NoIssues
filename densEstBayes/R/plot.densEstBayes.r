########## R function: plot.densEstBayes ##########

# For plotting a densEstBayes() fit object:

# Last changed: 14 SEP 2020

plot.densEstBayes <- function(x,plotIt=TRUE,credLev=0.95,gridSize=1001,varBand=TRUE,shade=TRUE,
                           estCol="darkgreen",varBandCol=NULL,axisCol="slateblue",
                           add=FALSE,lwd=2,xlab=NULL,ylab=NULL,...)
{
   fitObject <- x

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

   # Extract `xname':

   xname <- fitObject$xname

   # Set plotting matrices: 

   anchorLow <- -0.05  ;  anchorUpp <- 1.05
   xTrang <- seq(anchorLow,anchorUpp,length=gridSize)
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

   if (plotIt)
   {
      # Determine the default variability band colour if not specified:

      if (is.null(varBandCol))
      {
         if (!shade) varBandCol <- estCol
         if (shade) varBandCol <- "palegreen"
      } 

      # Determine the vertical range:

      if (!varBand)
         ylimVal <- range(densEstg)

      if (varBand)
         ylimVal <- range(c(densLowg,densUppg))

      # Do the plot:

      if (!add)
      {
         if (is.null(xlab)) xlab <- xname 
         if (is.null(ylab)) ylab <- "density"

         plot(0,type="n",bty="l",xlim=range(xg),
           ylim=ylimVal,xlab=xlab,ylab=ylab,...)
      } 

      if (varBand)
      {
         if (!shade)
         { 
            lines(xg,densLowg,col=varBandCol,lwd=lwd,lty=2)
            lines(xg,densUppg,col=varBandCol,lwd=lwd,lty=2)
         }
         if (shade)
            polygon(c(xg,rev(xg)),c(densLowg,rev(densUppg)),
                    col=varBandCol,border=FALSE)                  
      }
      lines(xg,densEstg,col=estCol,lwd=lwd)
      abline(0,0,col=axisCol)
      invisible()
   }
   if (!plotIt)
      return(list(xg=xg,densEstg=densEstg,densLowg=densLowg,densUppg=densUppg))
}

############ End of plot.densEstBayes ############
