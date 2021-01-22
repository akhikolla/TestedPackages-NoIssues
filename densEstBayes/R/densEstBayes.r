########## R function: densEstBayes ##########

# For performing density estimation via one of 
# four Bayesian inference engines.

# Last changed: 28 SEP 2020

densEstBayes <- function(x,method="slice",verbose=TRUE,
                         control=densEstBayes.control())
{
   # Extract `xname' for possible use in plot.densEstBayes():

   xname <- deparse(substitute(x))

   # Check the legality of x:

   x <- as.vector(x)
   if (length(x)<10) stop("x must have at least 10 entries.")

   # Check legality of method:

   if (!any(method==c("HMC","NUTS","slice","SMFVB")))
   {
      warnStr1 <- "The inputted method is not one of the available options."
      warnStr2 <- "The slice sampling default method was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      method <- "slice"
   }

   # Pass the verbose specification to the msgCode variable:

   msgCodeVerbose <- as.numeric(verbose)

   # Unpack the control values:

   range.x <- control$range.x
   numBins <- control$numBins
   numBasis <- control$numBasis
   sigmabeta <- control$sigmabeta
   ssigma <- control$ssigma
   convToler <- control$convToler
   maxIter <- control$maxIter
   nWarm <- control$nWarm
   nKept <- control$nKept
   nThin <- control$nThin
   msgCodeControl <- control$msgCode

   # Sort out conflicts between the "verbose" specification and the 
   # "msgCode" specification.

   if ((msgCodeVerbose==1)&(msgCodeControl==0))
   {              
      warnStr1 <- "The verbose and msgCode specifications conflict with"
      warnStr2 <- "each other. The default value of msgCode=1 was used."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      msgCodeControl <- 1 
   }

   # Determine the value of "msgCode" based on current values of 
   # "msgCodeVerbose" and "msgCodeControl".

   if (msgCodeVerbose==0) msgCode <- 0
   if (msgCodeVerbose>0) msgCode <- msgCodeControl
   
   # Set the default values of nWarm and nKept in the event
   # that no non-NULL values are specified:

   if (is.null(nWarm))
   {
      if (any(method==c("HMC","NUTS"))) nWarm <- 1000
      if (method=="slice") nWarm <- 100
   }
   if (is.null(nKept))
   {
      if (any(method==c("HMC","NUTS"))) nKept <- 1000
      if (method=="slice") nKept <- 1000
   }

   # Set lower and upper anchor values for spline basis 
   # functions after transformation to the unit interval:

   anchorLow <- -0.05  ;  anchorUpp <- 1.05

   # Set default value of range.x:

   if (is.null(range.x))
   {
      beyondDataFraction <- 0.05

      range.x <- c((1+beyondDataFraction)*min(x)-beyondDataFraction*max(x),
                   (1+beyondDataFraction)*max(x)-beyondDataFraction*min(x))
   } 

   # Check that range.x includes at least some of the data:

   if ((max(x)<=range.x[1])|(min(x)>=range.x[2]))
      stop("The inputted range.x vector is such that no data are contained\n
            within its interval.")

   # Discard data outside of the range.x interval:

   omitInds <- (1:length(x))[(x<range.x[1])|(x>range.x[2])]
   if (length(omitInds)>0)
      x <- x[-omitInds]

   # Obtain data transformed to the unit interval:

   xTran <- (x - range.x[1])/(range.x[2] - range.x[1])

   # Obtain Bayesian density estimate:
           
   gridPoints <- seq(anchorLow,anchorUpp,length=numBins)

   # Obtain grid counts:

   gridCounts <- round(linbin(xTran,gridPoints,TRUE))

   # Obtain knots and spline basis function design matrices:

   X <- cbind(1,gridPoints)
   intKnots <-  as.numeric(quantile(unique(xTran),seq(0,1,length=numBasis)[-c(1,numBasis)]))
   Z <- ZOSull(gridPoints,intKnots=intKnots,range.x=c(anchorLow,anchorUpp))
   ncZ <- ncol(Z)

   # Obtain density estimate:

   determFitObj <- NULL   ;   stochaFitObj <- NULL

   if (method=="SMFVB")
      determFitObj <- PoiSimMixModSMFVB(gridCounts,X,Z,sigmabeta,ssigma,
                                        convToler,maxIter)
 
   if (any(method==c("HMC","NUTS","slice")))
      stochaFitObj <- PoiSimMixModMCMC(gridCounts,X,Z,method,sigmabeta,
                                       ssigma,nWarm,nKept,nThin,msgCode)
 
   # Obtain sample hexiles for the possibly truncated and
   # transformed data (to be used for chain checks):

   sampHexTran <- quantile(xTran,(1:5)/6)  

   # Prepare and return the output object:

   outObj <- list(method=method,range.x=range.x,intKnots=intKnots,determFitObj=determFitObj,
                  stochaFitObj=stochaFitObj,xname=xname,sampHexTran=sampHexTran)

   class(outObj) <- "densEstBayes"

   return(outObj)
}   

############ End densEstBayes ###########

