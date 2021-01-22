/********** C++ function: linbin **********/

/* Performs linear binning of a set of real numbers
   over an equally-spaced grid. */

/* Last changed: 16 SEP 2020 */

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]

NumericVector linbin(arma::vec x,arma::vec gpoints,bool truncate)
{
   /* Declare all non-input variables: */

   int n = x.n_elem;
   int M = gpoints.n_elem;
   double delta;
   double lxi;
   int li;
   double rem;
   NumericVector gcounts(M);

   /* Set the grid width. */
  
   delta = (gpoints(M-1) - gpoints(0))/(M-1);
   
   for (int i = 0; i < n; i++)
   {
      lxi = ((x(i)-gpoints(0))/delta) + 1.0;
      li = floor(lxi);
      rem = lxi - li;

      if ((li>=1)&(li<M))
      {
         gcounts(li-1) =  gcounts(li-1) + 1.0 - rem;
         gcounts(li) = gcounts(li) + rem;
      }
      
      if ((li<1)&(not truncate))
         gcounts(0) = gcounts(0) + 1;

      if ((li>=M)&(not truncate))
         gcounts(M-1) = gcounts(M-1) + 1;
   }

   return gcounts;
}

/************ End of linbin ************/

