#define ARMA_DONT_PRINT_ERRORS
#include "func.h"
#include <RcppArmadillo.h>

using namespace Rcpp;

///////////////////////////////////
//         MEMBERSHIP DEGREE    //
//                              //
///////////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]

arma::mat memb_degree(arma::mat D, double m, unsigned int n, unsigned int k, unsigned int p)
{

  arma::mat out(n,k); out.zeros();
  arma::rowvec d(n); d.zeros();
  arma::uword mm = 0;

  for(int i=0; i<n;i++)
  {

    d = D.row(i);

    if(d.min() == 0)
    {

      mm = d.index_min();
      out(i,mm) = 1;

    }else{

      for(int j=0; j<k;j++)
      {

        out(i,j) = pow(1/arma::as_scalar(D(i,j)),1/(m-1.0)) / sum(pow(1/D.row(i),1/(m-1.0)));

      }
    }
  }

  return out;
}

/////////////////////////////////////////////////

///////////////////////////////////
//   MEMBERSHIP DEGREE   POL     //
///////////////////////////////////

// [[Rcpp::export]]

arma::mat memb_degree_pf(arma::mat D, double b, unsigned int n, unsigned int k, unsigned int p)
{

  arma::mat out(n,k); out.zeros();
  arma::rowvec d(n); d.zeros();
  unsigned int ki = 0;
  unsigned int kok = 0;
  arma::uword mm = 0;

  arma::uvec id;

  for(int i=0; i<n;i++)
  {

    d = D.row(i);

    if(d.min() == 0)
    {

      mm = d.index_min();
      out(i,mm) = 1;

    }else{

      ki = 0;
      kok = 0;
      d = sort(d);
      while(ki < k)

      {
        if(d(ki) * sum(1/d.subvec(0,ki)) <= ((1/b) + ki))
        {
          kok = ki;
          ki++;
        }else{

          ki++;
        }
      }
      out.row(i) = 1/(1-b)*((1+b*(kok))/(D.row(i)*sum(1/d.subvec(0,kok)))-b);
      out.row(i) = replace(out.row(i),0,0);

    }
  }

  return out;

}

/////////////////////////////////////////////////

///////////////////////////////////
//   MEMBERSHIP DEGREE ENTROPIC  //
///////////////////////////////////

// [[Rcpp::export]]

arma::mat memb_degree_ent(arma::mat D, double ent, unsigned int n, unsigned int k, unsigned int p)
{

  arma::mat out(n,k); out.zeros();
  arma::rowvec d(n); d.zeros();
  arma::uword mm = 0;

  double    eps = std::numeric_limits<double>::epsilon();

  for(int i=0; i<n;i++)
  {

    d = D.row(i);

    if(d.min() == 0)
    {

      mm = d.index_min();
      out(i,mm) = 1;

    }else{

      for(int j=0; j<k;j++)
      {
        out(i,j) = (exp(-arma::as_scalar(D(i,j))/ent)) / sum(exp(D.row(i)/(-ent)));
        if(arma::is_finite(out(i,j)) == false)
        {
          stop("Some membership degrees are NaN (Suggestion: run FKM.ent using standardized data)");
        }

        if(out(i,j) < eps)
        {
          out(i,j) = eps;
        }
      }
    }
  }

  return out;
}

/////////////////////////////////////////////////

///////////////////////////////////
//   MEMBERSHIP DEGREE MEDOID   //
///////////////////////////////////

// [[Rcpp::export]]

arma::mat memb_degree_medoid(arma::mat D, arma::uvec medoid, double m, unsigned int n, unsigned int k, unsigned int p)
{

  arma::mat out(n,k); out.zeros();
  arma::rowvec d(n); d.zeros();

  bool index_temp_med = true;
  arma::uvec find_temp(1); find_temp.zeros();

  for(int i=0; i<n;i++)
  {

    d = D.row(i);
    index_temp_med = Match(i,medoid);

    if(d.min() == 0)
    {

      arma::uword mm = d.index_min();
      out(i,mm) = 1;

    }else{

      for(int j=0; j<k;j++)
      {

        out(i,j) = pow(1/arma::as_scalar(D(i,j)),1/(m-1.0)) / sum(pow(1/D.row(i),1/(m-1.0)));
      }
    }
  }

  return out;
}

/////////////////////////////////////////////////
