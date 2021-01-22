#define ARMA_DONT_PRINT_ERRORS
#include "func.h"
#include <RcppArmadillo.h>

using namespace Rcpp;

///////////////////////////////////
//         CENTROIDS             //
///////////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
arma::mat centroids_FKM(arma::mat data, arma::mat U, unsigned int n, unsigned int k, unsigned int p, double m)
{

  arma::mat out(k,p, arma::fill::zeros);
  arma::mat U_new = trans(U);

  for(int i=0; i<k; i++)
  {

    out.row(i) = (pow(U_new.row(i),m) * data) / arma::as_scalar(sum(pow(U.col(i),m)));

  }

  return out;

}

/////////////////////////////////////////////////

// [[Rcpp::export]]

arma::mat centroids_FKM_ent(arma::mat data, arma::mat U, unsigned int n, unsigned int k, unsigned int p)
{

  arma::mat out(k,p, arma::fill::zeros);
  arma::mat U_new = trans(U);

  for(int i=0; i<k; i++)
  {

    out.row(i) = (U_new.row(i) * data) / arma::as_scalar(sum(U.col(i)));

  }

  return out;

}

/////////////////////////////////////////////////


// [[Rcpp::export]]
arma::mat centroids_FKM_pf(arma::mat data, arma::mat U, unsigned int n, unsigned int k, unsigned int p, double b)
{

  arma::mat out(k,p, arma::fill::zeros);
  arma::mat U_new = trans(U);

  double eps = std::numeric_limits<double>::epsilon();

  double bb1 = (1-b)/(1+b);
  double bb2 = (2*b) / (1+b);

  for(int i=0; i<k; i++)
  {

    out.row(i) = ((bb1*pow(U_new.row(i),2.0) + bb2 * U_new.row(i)) * data) /  (arma::as_scalar(sum(bb1 * pow(U.col(i),2.0) + bb2 * U.col(i))) + eps);

  }

  return out;

}

/////////////////////////////////////////////////


// [[Rcpp::export]]
List medoids_FKM(arma::mat data, arma::mat U, unsigned int n, unsigned int k, unsigned int p, double m)
{

  arma::mat out(k,p, arma::fill::zeros);
  double min_med_const = 1e5 * accu(pow(data,2));
  double min_med_old = 0;
  double min_med = 0;

  arma::uvec medoid(k); medoid.zeros();

  bool index_temp_med = true;
  bool value_temp_med = true;


  for(int c = 0; c < k; c++)
  {
    min_med_old = min_med_const;

    for(int i=0; i<n; i++)
    {
      min_med = 0;

      for(int j=0; j<n; j++)
      {
        min_med = min_med + pow(U(j,c),m) * sum(pow(data.row(i) - data.row(j),2.0));

      }

      index_temp_med = Match(i,medoid);
      value_temp_med = (min_med < min_med_old);

      if(index_temp_med && value_temp_med)
      {
        min_med_old = min_med;
        medoid(c) = i;
      }

    }

    out.row(c) = data.row(medoid(c));
  }
  return List::create(Rcpp::Named("H") = out,
                      Rcpp::Named("medoid") = medoid + 1);

}

/////////////////////////////////////////////////
