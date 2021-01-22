#define ARMA_DONT_PRINT_ERRORS
#include "func.h"
#include <RcppArmadillo.h>

using namespace Rcpp;



///////////////////////////////////
// EUCLIDEAN DISTANCE            //
///////////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

arma::mat euclidean_distance(arma::mat data, arma::mat H, unsigned int n, unsigned int k, bool Square = false)
{

  arma::mat out(n,k); out.zeros();

  for(int i=0; i<n; i++)
  {
    for(int j=0; j<k; j++)
    {

      out(i,j) = sum(square((data.row(i) - H.row(j))));
      if(Square == TRUE) {out(i,j) = sqrt(out(i,j));}

    }

  }


  return out;

}

/////////////////////////////////////////////////

// [[Rcpp::export]]

arma::mat euclidean_distance_gkb(arma::mat data, arma::mat H, arma::cube F, unsigned int n, unsigned int k, bool Square = false)
{

  arma::mat out(n,k); out.zeros();

  for(int i=0; i<n; i++)
  {
    for(int j=0; j<k; j++)
    {
      out(i,j) = arma::as_scalar(arma::dot((data.row(i)-H.row(j)) * arma::pinv(F.slice(j)) , (data.row(i)-H.row(j)).t()));

    }
  }

  return out;
}

/////////////////////////////////////////////////

// [[Rcpp::export]]

arma::mat euclidean_distance_gk(arma::mat data, arma::mat H, arma::cube F, arma::mat D_old, unsigned int n, unsigned int k, unsigned int p,bool Square = false)
{

  arma::mat out(n,k); out.zeros();
  arma::mat outReturn(n,k); outReturn.zeros();
  arma::mat inv_temp(p,p); inv_temp.zeros();
  arma::rowvec temp = data.row(1);
  bool check = false;

  for(int i=0; i<n; i++)
  {
    for(int j=0; j<k; j++)
    {
      inv_temp = InvCheck(F.slice(j));
      check = (inv_temp.is_empty() == true);
      if(check)
      {
        break;
      }

      out(i,j) = arma::as_scalar(arma::dot((data.row(i)-H.row(j)) * arma::pinv(F.slice(j)) , (data.row(i)-H.row(j)).t()));

    }
    if(check)
    {
      break;
    }
  }

  if(check)
  {
    outReturn = inv_temp;

  }else{
    outReturn = out;
  }
  return outReturn;
}

/////////////////////////////////////////////////


// [[Rcpp::export]]

arma::mat euclidean_distance_medoid(arma::mat data, arma::mat H, unsigned int n, unsigned int k, bool Square = false)
{
  arma::mat out(n,k); out.zeros();
  for(int i=0; i<n; i++)
  {
    for(int j=0; j<k; j++)
    {
      out(i,j) = sum(square((data.row(i) - H.row(j))));
      if(Square == true) {out(i,j) = sqrt(out(i,j));}
    }
  }
  return out;
}

/////////////////////////////////////////////////
