########## R function: PoiSimMixModSMFVB ##########

# For fitting and inference in a simple (single variance
# parameter) general design Poisson mixed model via 
# mean field variational Bayes.

# Last changed: 14 SEp 2020

PoiSimMixModSMFVB <- function(y,X,Z,sigmabeta,ssigma,convToler,maxIter)
{
   # Set required dimension variables:

   ncX <- ncol(X)  ;  ncZ <- ncol(Z)  ;  ncC <- ncX + ncZ

   # Set constant variational parameter values:

   kappa.q.sigsq <- 0.5*(ncZ + 1)
   kappa.q.a <- 1   

   # Set initial values of non-constant variational parameters:

   groupVec <- as.factor(rep(1,length(y)))
   fitPQL <- MASS::glmmPQL(y~-1+X,random=list(groupVec=nlme::pdIdent(~-1+Z)),family=poisson,
                              verbose=FALSE)

   mu.q.betau <- c(fitPQL$coef$fixed,unlist(fitPQL$coef$random))
   Sigma.q.betau <- diag(ncC)
   mu.q.recip.sigsq <- as.numeric(exp(-2*unlist(fitPQL$modelStruct)))

   # Form the C matrix:

   Cmat <- cbind(X,Z)

   # Perform iterations:

   converged <- FALSE
   itNum <- 0
   while (!converged)
   {
      itNum <- itNum + 1

      # Update q(beta,u) parameters:

      w.q.betau <- as.vector(exp(crossprod(t(Cmat),mu.q.betau) + 0.5*rowSums((crossprod(t(Cmat),Sigma.q.betau)*Cmat))))
      M.q.recip.sigsq <- diag(c(rep((1/sigmabeta^2),ncX),rep(mu.q.recip.sigsq,ncZ)))
      Sigma.q.betau <- try(solve(crossprod(Cmat,(w.q.betau*Cmat)) + M.q.recip.sigsq),silent=TRUE)

      # Check whether solve() successfuly or not:

      failedSolveCall <-  is.null(dim(Sigma.q.betau))
      if (failedSolveCall)
      { 
          errMsg1 <- "\n\n Call to densEstBayes() with method=\"SMFVB\" has failed due\n"
          errMsg2 <- "to a singular matrix and/or non-convergence problem. Use\n"
          errMsg3 <- "method=\"slice\", method=\"NUTS\" or method=\"HMC\" instead.\n\n"
          stop(paste(errMsg1,errMsg2,errMsg3),call.=FALSE) 
      } 

      if(!failedSolveCall)
      {      
         mu.q.betau <- mu.q.betau + crossprod(Sigma.q.betau,(crossprod(Cmat,(y-w.q.betau)) - M.q.recip.sigsq%*%mu.q.betau))

         # Update q(a) parameters:

         lambda.q.a <- mu.q.recip.sigsq + (1/ssigma^2)
         mu.q.recip.a <- kappa.q.a/lambda.q.a

         # Update q(sigsq) parameters:

         lambda.q.sigsq <- mu.q.recip.a + 0.5*(sum(mu.q.betau[-(1:ncX)]^2) + sum(diag(Sigma.q.betau)[-(1:ncX)]))

         mu.q.recip.sigsq <- kappa.q.sigsq/lambda.q.sigsq

         if (itNum>1)
            relErr <- abs((mu.q.recip.sigsq - mu.q.recip.sigsq.OLD)/mu.q.recip.sigsq.OLD)

         mu.q.recip.sigsq.OLD <- mu.q.recip.sigsq

         # Do convergence check:

         if (itNum>1)
           if (relErr<convToler)
               converged <- TRUE

         if (itNum>=maxIter)
         {
            converged <- TRUE
            warning("Semiparametric mean field variational Bayes did not converge after maximum number of iterations.\n")
         }
      }
   }

   # Return q-density parameters:

   return(list(mu.q.betau=mu.q.betau,Sigma.q.betau=Sigma.q.betau,
               kappa.q.sigsq=kappa.q.sigsq,lambda.q.sigsq=lambda.q.sigsq))
}  

############ End of PoiSimMixModSMFVB ############

