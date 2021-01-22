#define ARMA_DONT_PRINT_ERRORS
#include "func.h"
#include <RcppArmadillo.h>

using namespace Rcpp;



///////////////////////////////////
// COVARIANCE MATRICES            //
///////////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::cube F_gkb(arma::mat data, arma::mat U, arma::mat H, arma::mat F0, double m, double gam, unsigned int n, unsigned int k, int p, double mcn, arma::vec vp)
{

  //arma::vec dd(k); dd.zeros();
  double dd = 0;
  arma::cube F(p,p,k); F.zeros();
  double kappa = 0;
  bool kappaCond = true;
  bool kappaFin = true;

  arma::vec eigval;
  arma::mat eigvec;
  arma::mat eigvec1;

  arma::mat diag_temp(p,p); diag_temp.zeros();

  double em = 0;
  double em_mcm = 0;

  for(int i=0; i<k; i++)
  {

    for(int j=0; j<n; j++)
    {

      F.slice(i) += pow(U(j,i),m)* (data.row(j) - H.row(i)).t()*(data.row(j) - H.row(i));

    }

    F.slice(i) =  F.slice(i) / arma::as_scalar(sum(pow(U.col(i),m)));
    F.slice(i) = (1-gam) * F.slice(i) + gam * F0;

    kappa = arma::cond(F.slice(i));
    kappaCond = kappa > mcn;
    kappaFin = arma::is_finite(kappa);

    if(kappaCond | !kappaFin)
    {


      arma::svd(eigvec,eigval,eigvec1,F.slice(i), "std");
      em = max(eigval);
      em_mcm = em/mcn;
      eigval.elem(arma::find(eigval < em_mcm)).fill(em_mcm);
      diag_temp.diag() = eigval;
      F.slice(i) = eigvec * diag_temp * eigvec.i();

    }

    dd = std::abs(arma::det(F.slice(i)));
    F.slice(i) = F.slice(i) / (pow(dd, 1/(double)p) * vp(i));
  }

  return F;

}


/////////////////////////////////////////////////

// [[Rcpp::export]]
arma::cube F_gkb_ent(arma::mat data, arma::mat U, arma::mat H, arma::mat F0, double gam, unsigned int n, unsigned int k, int p, double mcn, arma::vec vp)
{

  double dd = 0;
  arma::cube F(p,p,k); F.zeros();
  double kappa = 0;
  bool kappaCond = true;
  bool kappaFin = true;

  arma::vec eigval;
  arma::mat eigvec;
  arma::mat eigvec1;

  arma::mat diag_temp(p,p); diag_temp.zeros();

  double em = 0;
  double em_mcm = 0;

  for(int i=0; i<k; i++)
  {

    for(int j=0; j<n; j++)
    {

      F.slice(i) += U(j,i)* (data.row(j) - H.row(i)).t()*(data.row(j) - H.row(i));

    }

    F.slice(i) =  F.slice(i) / arma::as_scalar(sum(U.col(i)));
    F.slice(i) = (1-gam) * F.slice(i) + gam * F0;

    kappa = arma::cond(F.slice(i));

    kappaCond = kappa > pow(10.0,15.0);
    kappaFin = arma::is_finite(kappa);


    if(kappaCond | !kappaFin)
    {

      arma::svd(eigvec,eigval,eigvec1,F.slice(i));
      em = max(eigval);
      em_mcm = em/mcn;
      eigval.elem(arma::find(eigval < em_mcm)).fill(em_mcm);
      diag_temp.diag() = eigval;
      F.slice(i) = eigvec * diag_temp * eigvec.i();

    }

    dd = std::abs(arma::det(F.slice(i)));
    F.slice(i) = F.slice(i) / (pow(dd, 1/(double)p) * vp(i));

  }

  return F;

}

/////////////////////////////////////////////////


// [[Rcpp::export]]

arma::cube F_gk(arma::mat data, arma::mat U, arma::mat H, double m,unsigned int n, unsigned int k, int p, arma::vec vp )
{
  double dd = 0;
  arma::cube F(p,p,k); F.zeros();

  for(int i = 0; i<k; i++)
  {
    for(int j=0; j<n; j++)
    {

      F.slice(i) += pow(U(j,i),m)* (data.row(j) - H.row(i)).t()*(data.row(j) - H.row(i));

    }

    F.slice(i) = F.slice(i) / arma::as_scalar(sum(pow(U.col(i),m)));
    dd = std::abs(det(F.slice(i)));
    F.slice(i) = F.slice(i) / (pow(dd, 1/(double)p) * vp(i));
  }
  return F;
}

/////////////////////////////////////////////////


// [[Rcpp::export]]

arma::cube F_gk_ent(arma::mat data, arma::mat U, arma::mat H,unsigned int n, unsigned int k, int p, arma::vec vp )
{
  double dd = 0;
  arma::cube F(p,p,k); F.zeros();

  for(int i = 0; i<k; i++)
  {
    for(int j=0; j<n; j++)
    {
      F.slice(i) += U(j,i)* (data.row(j) - H.row(i)).t()*(data.row(j) - H.row(i));
    }

    F.slice(i) =  F.slice(i) / arma::as_scalar(sum(U.col(i)));
    dd = std::abs(det(F.slice(i)));
    F.slice(i) = F.slice(i) / (pow(dd, 1/(double)p) * vp(i));
  }

  return F;
}

/////////////////////////////////////////////////
