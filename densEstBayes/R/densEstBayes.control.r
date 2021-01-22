########## R function: densEstBayes.control ##########

# Control function for densEstBayes().

# Last changed: 14 SEP 2020

densEstBayes.control <- function(range.x=NULL,numBins=401,numBasis=50,
                              sigmabeta=1e5,ssigma=1000,convToler=1e-5,
                              maxIter=500,nWarm=NULL,nKept=NULL,nThin=1,
                              msgCode=1)
{
   # Make sure that range.x is legal:

   if (!is.null(range.x))
   {
      if (length(range.x)!=2)
      {
         warnStr1 <- "The inputted range.x vector is not of length 2."
         warnStr2 <- "The default value was used instead."
         warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
         range.x <- NULL
      }
  
      if (range.x[1]>=range.x[2])
      {
         warnStr1 <- "The inputted range.x vector has first entry equal or exceeding"
         warnStr2 <- "the second entry. The default value was used instead."
         warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
         range.x <- NULL
      }
   }

   # Make sure that numBins is legal:

   numBins <- round(numBins)
   if (numBins<21)
   {
      warnStr1 <- "The inputted number of bins is negative or too low."
      warnStr2 <- "The default value of 401 was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      numBins <- 401
   }

   # Make sure that numBasis is legal:

   numBasis <- round(numBasis)
   if ((numBasis<5)|(numBasis>0.5*numBins))
   {
      warnStr1 <- "The inputted number of basis functions is illegal."
      warnStr2 <- "The default value of 50 was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      numBasis <- 50
   }

   # Make sure that sigmabeta is legal:
  
   if (sigmabeta<0)
   {
      warnStr1 <- "The inputted fixed effects scale hyperparameter is negative."
      warnStr2 <- "The default value of 100000 was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      sigmabeta <- 100000
   }

   # Make sure that ssigma is legal:
  
   if (ssigma<0)
   {
      warnStr1 <- "The inputted standard deviation scale hyperparameter is negative."
      warnStr2 <- "The default value of 1000 was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      ssigma <- 1000
   }

   # Make sure that convToler is legal:

   if ((convToler<0)|(convToler>0.1))
   {
      warnStr1 <- "The inputted semiparametric mean field variational"
      warnStr2 <- "Bayes tolerance value is negative or exceeds 0.1."
      warnStr3 <- "The default value of 0.00001 was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n  ",warnStr3,"\n",sep=""),
              immediate.=TRUE)
      convToler <- 1e-5
   }

   # Make sure that maxIter is legal:

   maxIter <- round(maxIter)
   if (maxIter<10)
   {
      warnStr1 <- "The inputted maximum number of semiparametric mean field"
      warnStr2 <- "variational Bayes iterations is negative or too low."
      warnStr3 <- "The default value of 500 was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n  ",warnStr3,"\n",sep=""),
              immediate.=TRUE)
      maxIter <- 500
   }
   
   if (!is.null(nWarm))
   {
      # Make sure that nWarm is legal:

      nWarm <- round(nWarm)
      if (nWarm<0) 
      {
         warnStr1 <- "The inputted number of warmup Markov chain Monte Carlo iterations is negative."
         warnStr2 <- "The default value of 100 was used instead."
         warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
         nWarm <- 100
      }
   }

   # Make sure that nKept is legal:

   if (!is.null(nKept))
   {
      nKept <- round(nKept)
      if (nKept<0) 
      {
         warnStr1 <- "The inputted number of kept Markov chain Monte Carlo iterations is negative."
         warnStr2 <- "The default value of 1000 was used instead."
         warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
         nKept <- 1000
      }
   }

   # Make sure that nThin is legal:

   nThin <- round(nThin)
   if (nThin<0) 
   {
      warnStr1 <- "The inputted Markov chain Monte Carlo thinning factor is negative."
      warnStr2 <- "The default value of 1 was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      nThin <- 1
   }
   
   if (!is.null(nKept))
   {
      if (nThin>nKept) 
      {
         warnStr1 <- "The inputted Markov chain Monte Carlo thinning factor exceeds the"
         warnStr2 <- "number of kept values. The default value of 1 was used instead."
         warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
         nThin <- 1
      }
   }

   # Make sure that msgCode is legal:

   msgCode <- round(msgCode)
   if (!any(msgCode==(0:2)))
   {
      warnStr1 <- "The inputted message code number is not 0, 1 or 2."
      warnStr2 <- "The default value of 1 was used instead."
      warning(paste(warnStr1,"\n  ",warnStr2,"\n",sep=""),immediate.=TRUE)
      msgCode <- 1
   }

   return(list(range.x=range.x,numBins=numBins,numBasis=numBasis,
               sigmabeta=sigmabeta,ssigma=ssigma,convToler=convToler,
               maxIter=maxIter,nWarm=nWarm,nKept=nKept,nThin=nThin,
               msgCode=msgCode))
}

############ End of densEstBayes.control ############
