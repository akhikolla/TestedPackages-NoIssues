/********** C++ function: logUnnDens **********/

/* Computes the logarithm of the unnormalised density
   function that arises in slice sampling for the simple
   Poisson mixed model. */

/* Last changed: 17 SEP 2020 */

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]

double logUnnDens(int j,arma::vec betauj,arma::vec betauNotj,
                  arma::mat Cmat,arma::vec CTy,arma::vec sigsqbetau)
{
   /* Declare all non-input variables: */

   int ncC = Cmat.n_cols;
   arma::vec betau;
   arma::vec sttPc;
   arma::vec endPc;
   double firstTerm;
   double secondTerm;
   double thirdTerm;
   double theAnswer;

   /* Compute the first term: */

   firstTerm = CTy(j-1)*betauj(0);

   /* Compute the second term: */

   if (j==1) {betau = arma::join_cols(betauj,betauNotj);};
   if ((j>1)&(j<ncC))
   {
      sttPc = betauNotj.rows(0,j-2);
      endPc = arma::join_cols(betauj,betauNotj.rows(j-1,ncC-2));
      betau = arma::join_cols(sttPc,endPc);
   };
   if (j==ncC) {betau = arma::join_cols(betauNotj,betauj);};  
   secondTerm = -sum(exp(Cmat * betau));

   /* Compute the third term: */
 
   thirdTerm = -0.5*betauj(0)*betauj(0)/sigsqbetau(j-1);
 
   /* Combine terms to obtain the answer: */

   theAnswer = firstTerm + secondTerm + thirdTerm;

   return theAnswer;
}

/************ End of logUnnDens ************/

