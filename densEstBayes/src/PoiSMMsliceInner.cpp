/********** C++ function: PoiSMMsliceInner **********/

/* Carries out Markov chain Monte Carlo for a Normal  
   random sample.                                    */

/* Last changed: 17 SEP 2020 */

#include <RcppArmadillo.h>
#include "printPercMsgs.h"
#include "logUnnDens.h"

using namespace Rcpp;

// [[Rcpp::export]]

List PoiSMMsliceInner(int numMCMC,int ncX,arma::vec y,arma::mat Cmat,
                      arma::vec CTy,double sigmabetaHYP,double AHYP,          
                      int msgCode)
{
   /* Declare all non-input variables: */
                       
   int ncC = Cmat.n_cols;
   arma::mat betau(ncC,numMCMC);
   arma::vec sigsq(numMCMC);
   arma::vec sigsqbetau(ncC);
   double recipSqAHYP;
   int percCnt; 
   bool newPointFound;
   arma::vec betauCurr;
   arma::vec betauCurrj;
   arma::vec betauCurrNotj;
   double logUnnDensCurr;
   double aOnly; 
   double w;
   double unifVal;
   arma::vec betauNewj;
   arma::vec Lbetauj;
   arma::vec Ubetauj;
   double currThresh;
   double sumDists; 
   arma::vec LbetaujTilde;
   arma::vec UbetaujTilde;
   double shapeVal;
   double scaleVal;
   double aux;
   double unorm;

   /* Do some dummy initialisations to ensure
      correct allocation of space: */

   betauNewj = y.rows(1,1);
   Lbetauj = y.rows(1,1);
   Ubetauj = y.rows(1,1);
   LbetaujTilde = y.rows(1,1);
   UbetaujTilde = y.rows(1,1);

   /* Set constant quantity: */

   recipSqAHYP = 1.0/(AHYP*AHYP);

   /* Initialise chains: */

   sigsqbetau = arma::ones(ncC);
   for (int k = 0; k < ncX; k++)
   {
      sigsqbetau(k) = sigmabetaHYP*sigmabetaHYP;
   }

   betau.col(0) = 0.0*arma::ones(ncC);
   sigsq(0) = 1.0;

   /* Perform Markov chain Monte Carlo sampling: */

   sumDists = 0.0;
   percCnt = 0; 
   if (msgCode>0)
   {
      Rcout << "   The percentage of slice sampling completed is:" <<"\n";
      Rcout << "   ";
   };
   for (int gPos = 2; gPos <=numMCMC; gPos++)
   {
      /* Print percentage progess message. */

      percCnt = printPercMsgs(msgCode,numMCMC,gPos,percCnt);

      /* Draw sample for (beta,u) using slice sampling: */

      betauCurr = betau.col(gPos-2);
      
      for (int j = 1; j <= ncC; j++)
      {  
 
         betauCurrj = betauCurr.rows((j-1),(j-1)); 
        
         if (j==1){betauCurrNotj = betauCurr.rows(1,ncC-1);}; 
         if ((j>1)&(j<ncC))
         {
            betauCurrNotj = arma::join_cols(betauCurr.rows(0,j-2),betauCurr.rows(j,ncC-1));
         }; 
         if (j==ncC){betauCurrNotj = betauCurr.rows(0,ncC-2);}; 


         logUnnDensCurr = logUnnDens(j,betauCurrj,betauCurrNotj,Cmat,CTy,sigsqbetau);

         aOnly = logUnnDensCurr + log(Rcpp::runif(1)(0));  

         if (gPos==2) {w = 1.0;};
         if (gPos>2) {w = sumDists/(gPos-2);};

         Lbetauj(0) = betauCurrj(0) - w*Rcpp::runif(1)(0);    
         Ubetauj(0) = Lbetauj(0) + w;

         currThresh = logUnnDens(j,Lbetauj,betauCurrNotj,Cmat,CTy,sigsqbetau);
         while (aOnly<currThresh)
         {
            Lbetauj(0) = Lbetauj(0) - w; 
            currThresh = logUnnDens(j,Lbetauj,betauCurrNotj,Cmat,CTy,sigsqbetau);
         }

         currThresh = logUnnDens(j,Ubetauj,betauCurrNotj,Cmat,CTy,sigsqbetau);
         while (aOnly<currThresh)
         {
            Ubetauj(0) = Ubetauj(0) + w;
            currThresh = logUnnDens(j,Ubetauj,betauCurrNotj,Cmat,CTy,sigsqbetau);
         }

         LbetaujTilde(0) = Lbetauj(0);
         UbetaujTilde(0) = Ubetauj(0);
         newPointFound = FALSE;         
         while (not newPointFound)
         {  
            unifVal = Rcpp::runif(1)(0);  
            betauNewj(0) = LbetaujTilde(0) + unifVal*(UbetaujTilde(0) - LbetaujTilde(0)); 
            if (aOnly<logUnnDens(j,betauNewj,betauCurrNotj,Cmat,CTy,sigsqbetau))
            {
               betau(j-1,gPos-1) = betauNewj(0);
               betauCurr(j-1) = betauNewj(0);
               newPointFound = TRUE;
            };             
            if (betauNewj(0)<betauCurrj(0))  {LbetaujTilde(0) = betauNewj(0);};
            if (betauNewj(0)>=betauCurrj(0)) {UbetaujTilde(0) = betauNewj(0);};
         }   
         sumDists = sumDists + std::abs(betauNewj(0)-betauCurrj(0));   
      }
 
      /* Draw sample for the reciprocal of sigma squared: */

      shapeVal = 1.0; scaleVal = 1.0;
      aux = ((1.0/sigsq(gPos-2))+ recipSqAHYP)/(R::rgamma(shapeVal,scaleVal));

      unorm = 0.0;
      for (int k = ncX; k < ncC; k++)
      {
         unorm = unorm + betau(k,gPos-1)*betau(k,gPos-1);
      }  
      shapeVal = 0.5*(ncC - ncX + 1.0);
      
      sigsq(gPos-1) = ((1.0/aux) + 0.5*unorm)/(R::rgamma(shapeVal,scaleVal));  

      for (int k = ncX; k < ncC; k++)
      {
         sigsqbetau(k) = sigsq(gPos-1);
      }  
   }
   if (msgCode>0) {Rcout << "\n";};

   List MCMCsamples = List::create(
                           Named("betau",betau),
                           Named("sigsq",sigsq));

   return MCMCsamples;
}

/************ End of PoiSMMsliceInner ************/

