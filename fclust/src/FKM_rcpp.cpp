#define ARMA_DONT_PRINT_ERRORS
#include "func.h"
#include "centroids.h"
#include "memDeg.h"
#include "distance.h"
#include "F.h"

#include <RcppArmadillo.h>



// using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


//////////////////////////////////////////////
///         FUZZY ALGORITHMS                //
//     core function exported in R          //
//       to run fuzzy algorithms            //
//////////////////////////////////////////////

///////////////////////////////////
//              FKM              //
///////////////////////////////////


// [[Rcpp::export]]

List mainFKM(arma::mat data,
                    double m,
                    unsigned int n,
                    unsigned int p,
                    unsigned int k,
                    unsigned int rs,
                    double conv,
                    unsigned int maxit,
                    std::string index,
                    double alpha) {

  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U(n,k); U.zeros();
  arma::mat U_old(n,k); U_old.zeros();
  double func = 0;
  double func_opt = 0;

  double ind = 0;
  double ind_max = 0;

  arma::mat H_opt(k, p); H.zeros();
  arma::mat U_opt(n,k); U.zeros();

  for(int r=0; r<rs; r++)
  {



    int iter = 0;
    U = unifInit(n,k);
    U_old = U;

    bool prova = true;

    while(prova && (iter < maxit))
    {

      iter++;
      U_old = U;

      H = centroids_FKM(data, U, n, k, p, m);

      D = euclidean_distance(data, H, n, k);

      U = memb_degree(D, m, n, k, p);

      prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

    }

    func = accu(pow(U,m)%D);
    it(r) = iter;
    value(r) = func;
    if(arma::is_finite(func) == true)
    {
      if ( (r == 0) | (func < func_opt) )
      {
        U_opt=U;
        H_opt=H;
        func_opt=func;
      }
    }
  }

  ind = indices(index, data, U_opt, H_opt, m,  n, k, p,  exp(1.0), alpha);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

    if(index == "PE" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }

  return List::create(Rcpp::Named("H") = H_opt,
                      Rcpp::Named("U") = U_opt,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k);

}

// [[Rcpp::export]]

List mainFKM_U(arma::mat data,
              double m,
              unsigned int n,
              unsigned int p,
              unsigned int k,
              arma::mat U,
              double conv,
              unsigned int maxit,
              std::string index,
              double alpha) {

  int iter = 0;
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U_old = U;

  bool prova = true;
  double ind = 0;
  double ind_max = 0;
  double value = 0;

  while(prova && (iter < maxit))
  {

    iter++;
    U_old = U;

    H = centroids_FKM(data, U, n, k, p, m);

    D = euclidean_distance(data, H, n, k);

    U = memb_degree(D, m, n, k, p);

    prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

  }

  value = accu(pow(U,m)%D);

  ind = indices(index, data, U, H, m,  n, k, p,  exp(1.0), alpha);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }
  return List::create(Rcpp::Named("H") = H,
                      Rcpp::Named("U") = U,
                      Rcpp::Named("iter") = iter,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k);

}

//////////////////////////////////////////////////////////////////////////




///////////////////////////////////
//     ENTROPIC FKM              //
///////////////////////////////////





// [[Rcpp::export]]

List mainFKM_ent(arma::mat data,
                 double ent,
                 unsigned int n,
                 unsigned int p,
                 unsigned int k,
                 unsigned rs,
                 double conv,
                 unsigned int maxit,
                 std::string index,
                 double alpha) {

  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U(n,k); U.zeros();
  arma::mat U_old(n,k); U_old.zeros();

  double func = 0;
  double func_opt = 0;

  double ind = 0;
  double ind_max = 0;

  arma::mat H_opt(k, p); H.zeros();
  arma::mat U_opt(n,k); U.zeros();

  for(int r=0; r<rs; r++)
  {



    int iter = 0;
    U = unifInit(n,k);
    U_old = U;


    bool prova = true;

    while(prova && (iter < maxit))
    {

      iter++;
      U_old = U;

      H = centroids_FKM_ent(data, U, n, k, p);

      D = euclidean_distance(data, H, n, k);

      U = memb_degree_ent(D, ent, n, k, p);

      prova = arma::as_scalar(accu(abs(U_old - U))) > conv;

    }

    func = accu(U%D) + ent * accu(U % log(U));

    it(r) = iter;
    value(r) = func;

    if(arma::is_finite(func) == true)
    {
      if ((r == 0) | (func < func_opt))
      {
        U_opt=U;
        H_opt=H;
        func_opt=func;
      }
    }


  }


  ind = indices(index, data, U_opt, H_opt, 2,  n, k, p,  exp(1.0), alpha);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }

  return List::create(Rcpp::Named("H") = H_opt,
                      Rcpp::Named("U") = U_opt,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k);
}


// [[Rcpp::export]]

List mainFKM_ent_U(arma::mat data,
                   double ent,
                   unsigned int n,
                   unsigned int p,
                   unsigned int k,
                   arma::mat U,
                   double conv,
                   unsigned int maxit,
                   std::string index,
                   double alpha) {

  int iter = 0;
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U_old = U;

  bool prova = true;

  double value = 0;
  double ind = 0;
  double ind_max = 0;

  while(prova && (iter < maxit))
  {

    iter++;
    U_old = U;

    H = centroids_FKM_ent(data, U, n, k, p);

    D = euclidean_distance(data, H, n, k);

    U = memb_degree_ent(D, ent, n, k, p);

    prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

  }

  value = accu(U%D) + ent * accu(U % log(U));

  ind = indices(index, data, U, H, 2,  n, k, p,  exp(1.0), alpha);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }
  return List::create(Rcpp::Named("H") = H,
                      Rcpp::Named("U") = U,
                      Rcpp::Named("iter") = iter,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k);

}


///////////////////////////////////
//             NOISE FKM        //
///////////////////////////////////


// [[Rcpp::export]]

List mainFKM_noise(arma::mat data,
                   double m,
                   double delta,
                   unsigned int n,
                   unsigned int p,
                   unsigned int k,
                   unsigned rs,
                   double conv,
                   unsigned int maxit,
                   std::string index,
                   double alpha) {

  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U(n,k); U.zeros();
  arma::mat U_old(n,k); U_old.zeros();

  double func = 0;
  double func_opt = 0;
  double ind = 0;
  double ind_max = 0;

  arma::mat H_opt(k, p); H.zeros();
  arma::mat U_opt(n,k); U.zeros();

  int nan_check = 1;

  for(int r=0; r<rs; r++)
  {



    int iter = 0;
    U = unifInit(n,k);
    U_old = U;

    arma::rowvec d(n); d.zeros();
    arma::vec Uout(n); Uout.zeros();
    double q1 = 0;
    double q2 = 0;

    double deltasq = pow(delta,2.0);

    bool prova = true;

    while(prova && (iter < maxit))
    {

      iter++;
      U_old = U;

      H = centroids_FKM(data, U, n, k, p, m);

      D = euclidean_distance(data, H, n, k);

      for(int i=0; i<n;i++)
      {
        nan_check = arma::is_finite(U.row(i));
        if(nan_check == 0){
          break;
        }
        d = D.row(i);

        if(d.min() == 0)
        {

          arma::uword mm = d.index_min();
          U.row(i).zeros();
          U(i,mm) = 1;

        }else{

          for(int j=0; j<k;j++)
          {

            q1 = pow(1/arma::as_scalar(D(i,j)),1/(m-1.0)) / sum(pow(1/D.row(i),1/(m-1.0)));
            q2 = pow(arma::as_scalar(D(i,j)),1/(m-1.0)) / pow(delta,2/(m-1.0));
            U(i,j) = 1/((1/q1)+q2);

          }
        }

        Uout(i) = 1 - sum(U.row(i));
        Uout(i) = pow(Uout(i),m) * deltasq;

      }

      prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

    }

    if(nan_check == 0)
    {
      func = arma::datum::nan;
    }else{

      func = accu(pow(U,m)%D) + sum(Uout);
    }

    it(r) = iter;
    value(r) = func;

    if(arma::is_finite(func) == true)
    {
      if ((r == 0) | (func < func_opt))
      {
        U_opt=U;
        H_opt=H;
        func_opt=func;
      }
    }
  }


  if(arma::is_finite(func_opt)){

    ind = indices(index, data, U_opt, H_opt, m,  n, k, p,  exp(1.0), alpha);

    // if(index == "PE" || index == "PC" || index == "XB")
    //   {ind_max = -ind;
    //   }else{
    //     ind_max = ind;
    //   }

    if(index == "PE" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }
  }else{

    ind = arma::datum::nan;
    warning("Some membership degrees are NaN (Suggestion: Increase the number of starting points RS)");

  }

  return List::create(Rcpp::Named("H") = H_opt,
                      Rcpp::Named("U") = U_opt,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k);

}

// [[Rcpp::export]]

List mainFKM_noise_U(arma::mat data,
                     double m,
                     double delta,
                     unsigned int n,
                     unsigned int p,
                     unsigned int k,
                     arma::mat U,
                     double conv,
                     unsigned int maxit,
                     std::string index,
                     double alpha) {

  int iter = 0;
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U_old = U;

  arma::rowvec d(n); d.zeros();
  arma::vec Uout(n); Uout.zeros();
  arma::uword mm = 0;

  double q1 = 0;
  double q2 = 0;

  double ind = 0;
  double ind_max = 0;

  double deltasq = pow(delta,2.0);

  bool prova = true;

  double value = 0;

  int nan_check = 1;


  while(prova && (iter < maxit))
  {

    iter++;
    U_old = U;

    H = centroids_FKM(data, U, n, k, p, m);

    D = euclidean_distance(data, H, n, k);

    for(int i=0; i<n;i++)
    {

      d = D.row(i);

      if(d.min() == 0)
      {

        mm = d.index_min();
        U.row(i).zeros();
        U(i,mm) = 1;

      }else{

        for(int j=0; j<k;j++)
        {

          q1 = pow(1/arma::as_scalar(D(i,j)),1/(m-1.0)) / sum(pow(1/D.row(i),1/(m-1.0)));
          q2 = pow(arma::as_scalar(D(i,j)),1/(m-1.0)) / pow(delta,2/(m-1.0));
          nan_check = arma::is_finite(U.row(i));
          if(nan_check == 0){
            break;
          }
          U(i,j) = 1/((1/q1)+q2);

        }
      }

      Uout(i) = 1 - sum(U.row(i));
      Uout(i) = pow(Uout(i),m) * deltasq;

    }

    prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

  }

  if(nan_check == 0)
  {
    ind = arma::datum::nan;
    warning("Some membership degrees are NaN (Suggestion: Increase the number of starting points RS)");


  }else{

    value = accu(pow(U,m)%D) + sum(Uout);
    ind = indices(index, data, U, H, m,  n, k, p,  exp(1.0), alpha);

    // if(index == "PE" || index == "PC" || index == "XB")
    //   {ind_max = -ind;
    //   }else{
    //     ind_max = ind;
    //   }

    if(index == "PE" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }
  }

  return List::create(Rcpp::Named("H") = H,
                      Rcpp::Named("U") = U,
                      Rcpp::Named("iter") = iter,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k);

}

//////////////////////////////////////////////////////////////////////////

///////////////////////////////////
//      POLINOMIAL FKM           //
///////////////////////////////////







// [[Rcpp::export]]

List mainFKM_pf(arma::mat data,
                double b,
                unsigned int n,
                unsigned int p,
                unsigned int k,
                unsigned rs,
                double conv,
                unsigned int maxit,
                std::string index,
                double alpha) {

  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U(n,k); U.zeros();
  arma::mat U_old(n,k); U_old.zeros();

  double func = 0;
  double func_opt = 0;
  double ind = 0;
  double ind_max = 0;

  arma::mat H_opt(k, p); H.zeros();
  arma::mat U_opt(n,k); U.zeros();

  for(int r=0; r<rs; r++)
  {



    int iter = 0;
    U = unifInit(n,k);
    U_old = U;

    bool prova = true;

    while(prova && (iter < maxit))
    {

      iter++;
      U_old = U;

      H = centroids_FKM_pf(data, U, n, k, p, b);

      D = euclidean_distance(data, H, n, k);

      U = memb_degree_pf(D, b, n, k, p);

      prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

    }

    func = accu((((1-b)/(1+b))*pow(U,2.0) + ((2*b)/(1+b))*U) % D);
    it(r) = iter;
    value(r) = func;

    if(arma::is_finite(func) == true)
    {
      if ((r == 0) | (func < func_opt))
      {
        U_opt=U;
        H_opt=H;
        func_opt=func;
      }
    }



  }


  ind = indices(index, data, U_opt, H_opt, 2,  n, k, p,  exp(1.0), alpha);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }

  return List::create(Rcpp::Named("H") = H_opt,
                      Rcpp::Named("U") = U_opt,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k);

}

// [[Rcpp::export]]

List mainFKM_pf_U(arma::mat data,
                  double b,
                  unsigned int n,
                  unsigned int p,
                  unsigned int k,
                  arma::mat U,
                  double conv,
                  unsigned int maxit,
                  std::string index,
                  double alpha) {

  int iter = 0;
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U_old = U;

  bool prova = true;
  double ind = 0;
  double ind_max = 0;
  double value = 0;

  while(prova && (iter < maxit))
  {

    iter++;
    U_old = U;

    H = centroids_FKM_pf(data, U, n, k, p, b);

    D = euclidean_distance(data, H, n, k);

    U = memb_degree_pf(D, b, n, k, p);

    prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

  }

  value = accu((((1-b)/(1+b))*pow(U,2.0) + ((2*b)/(1+b))*U) % D);

  ind = indices(index, data, U, H, 2,  n, k, p,  exp(1.0), alpha);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }

  return List::create(Rcpp::Named("H") = H,
                      Rcpp::Named("U") = U,
                      Rcpp::Named("iter") = iter,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k);

}


//////////////////////////////////////////////////////////////////////

///////////////////////////////////
//   ENTROPIC - NOISE FKM       //
///////////////////////////////////

// [[Rcpp::export]]

List mainFKM_ent_noise(arma::mat data,
                       double ent,
                       double delta,
                       unsigned int n,
                       unsigned int p,
                       unsigned int k,
                       unsigned rs,
                       double conv,
                       unsigned int maxit,
                       std::string index,
                       double alpha) {

  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U(n,k); U.zeros();
  arma::mat U_old(n,k); U_old.zeros();

  arma::rowvec d(n); d.zeros();
  arma::uword mm = 0;

  double func = 0;
  double func_opt = 0;
  double ind = 0;
  double ind_max = 0;

  arma::mat H_opt(k, p); H.zeros();
  arma::mat U_opt(n,k); U.zeros();

  arma::vec Uout(n); Uout.zeros();

  double eps = std::numeric_limits<double>::epsilon();

  for(int r=0; r<rs; r++)
  {



    int iter = 0;
    U = unifInit(n,k);
    U_old = U;

    bool prova = true;

    while(prova && (iter < maxit))
    {

      iter++;
      U_old = U;

      H = centroids_FKM_ent(data, U, n, k, p);

      D = euclidean_distance(data, H, n, k);

      for(int i=0; i<n;i++)
      {

        d = D.row(i);

        if(d.min() == 0)
        {

          mm = d.index_min();
          U.row(i).zeros();
          U(i,mm) = 1;

        }else{

          for(int j=0; j<k;j++)
          {
            U(i,j) = exp(-arma::as_scalar(D(i,j))/ent) / (sum(exp(D.row(i)/(-ent))) + exp(- pow(delta,2.0) / ent));

            if(arma::is_finite(U(i,j)) == false)
            {
              stop("Some membership degrees are NaN (Suggestion: run FKM.ent.noise using standardized data)");
            }

            if(U(i,j) < eps)
            {
              U(i,j) = eps;
            }
          }
          Uout(i) = 1 - sum(U.row(i));
        }
      }

      prova = arma::as_scalar(accu(abs(U_old - U))) > conv;

    }

    func = accu(U%D) + ent * accu(U % log(U)) + sum(Uout) * pow(delta,2.0);

    it(r) = iter;
    value(r) = func;

    if(arma::is_finite(func) == true)
    {
      if ((r == 0) | (func < func_opt))
      {
        U_opt=U;
        H_opt=H;
        func_opt=func;
      }
    }


  }


  ind = indices(index, data, U_opt, H_opt, 2,  n, k, p,  exp(1.0), alpha);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }

  return List::create(Rcpp::Named("H") = H_opt,
                      Rcpp::Named("U") = U_opt,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k);

}


// [[Rcpp::export]]

List mainFKM_ent_noise_U(arma::mat data,
                         double ent,
                         double delta,
                         unsigned int n,
                         unsigned int p,
                         unsigned int k,
                         arma::mat U,
                         double conv,
                         unsigned int maxit,
                         std::string index,
                         double alpha) {

  int iter = 0;
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U_old = U;
  arma::uword mm = 0;

  arma::rowvec d(n); d.zeros();

  bool prova = true;
  double ind = 0;
  double ind_max = 0;

  double value = 0;

  arma::vec Uout(n); Uout.zeros();

  double eps = std::numeric_limits<double>::epsilon();

  while(prova && (iter < maxit))
  {

    iter++;
    U_old = U;

    H = centroids_FKM_ent(data, U, n, k, p);

    D = euclidean_distance(data, H, n, k);

    for(int i=0; i<n;i++)
    {

      d = D.row(i);

      if(d.min() == 0)
      {

        mm = d.index_min();
        U.row(i).zeros();
        U(i,mm) = 1;

      }else{

        for(int j=0; j<k;j++)
        {
          U(i,j) = exp(-arma::as_scalar(D(i,j))/ent) / (sum(exp(D.row(i)/(-ent))) + exp(- pow(delta,2.0) / ent));

          if(arma::is_finite(U(i,j)) == false)
          {
            stop("Some membership degrees are NaN (Suggestion: run FKM.ent.noise using standardized data)");
          }

          if(U(i,j) < eps)
          {
            U(i,j) = eps;
          }
        }
        Uout(i) = 1 - sum(U.row(i));
      }
    }

    prova = arma::as_scalar(accu(abs(U_old - U))) > conv;

  }



  value = accu(U%D) + ent * accu(U % log(U)) + sum(Uout) * pow(delta,2.0);


  ind = indices(index, data, U, H, 2,  n, k, p,  exp(1.0), alpha);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }

  return List::create(Rcpp::Named("H") = H,
                      Rcpp::Named("U") = U,
                      Rcpp::Named("iter") = iter,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k);

}

///////////////////////

///////////////////////////////////
//POLYNOMIAL - NOISE FKM         //
///////////////////////////////////

// [[Rcpp::export]]

List mainFKM_pf_noise(arma::mat data,
                      double b,
                      double delta,
                      unsigned int n,
                      unsigned int p,
                      unsigned int k,
                      unsigned rs,
                      double conv,
                      unsigned int maxit,
                      std::string index,
                      double alpha) {

  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U(n,k); U.zeros();
  arma::mat U_old(n,k); U_old.zeros();

  double func = 0;
  double func_opt = 0;
  double ind = 0;
  double ind_max = 0;

  arma::mat H_opt(k, p); H.zeros();
  arma::mat U_opt(n,k); U.zeros();

  arma::rowvec d(n); d.zeros();
  arma::uword mm = 0;
  unsigned int ki = 0;
  unsigned int kok = 0;
  arma::vec Uout(n); Uout.zeros();
  arma::uvec id;

  arma::rowvec deltatemp(1); deltatemp = delta;

  int nan_check = 1;


  for(int r=0; r<rs; r++)
  {



    int iter = 0;
    U = unifInit(n,k);
    U_old = U;

    bool prova = true;

    double deltasq = pow(delta,2.0);

    while(prova && (iter < maxit))
    {

      iter++;
      U_old = U;

      H = centroids_FKM_pf(data, U, n, k, p, b);

      D = euclidean_distance(data, H, n, k);

      for(int i=0; i<n;i++)
      {
        nan_check = arma::is_finite(U.row(i));
        if(nan_check == 0){
          break;
        }
        d = D.row(i);

        if(d.min() == 0)
        {

          mm = d.index_min();
          U.row(i).zeros();
          U(i,mm) = 1;

        }else{

          ki = 0;
          kok = 0;
          d = arma::join_rows(d,deltatemp);
          d = sort(d);
          while(ki < k)

          {
            //if((d(ki) * (sum(1/d.subvec(0,ki)) + (1/deltasq))) <= ((1/b) + ki))
            if((d(ki) * (sum(1/d.subvec(0,ki)))) <= ((1/b) + ki))
            {
              kok = ki;
              ki++;
            }else{

              ki++;
            }
          }
          U.row(i) = 1/(1-b)*(((1+b*(kok))/(D.row(i)*(sum(1/d.subvec(0,kok)))))-b);;
          U.row(i) = replace(U.row(i),0,0);

        }

        Uout(i) = 1 - sum(U.row(i));
      }

      prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

    }


    if(nan_check == 0)
    {
      func = arma::datum::nan;
    }else{

      func = accu((((1-b)/(1+b))*pow(U,2.0) + ((2*b)/(1+b))*U) % D) + accu((((1-b)/(1+b))*pow(U,2.0)* deltasq)) + sum(((2*b)/(1+b))*Uout * deltasq);
    }

    it(r) = iter;
    value(r) = func;

    if(arma::is_finite(func) == true)
    {
      if ((r == 0) | (func < func_opt))
      {
        U_opt=U;
        H_opt=H;
        func_opt=func;
      }
    }
  }

  ind = 1;
  ind_max = ind;
  if(arma::is_finite(func_opt)){

    ind = indices(index, data, U_opt, H_opt, 1,  n, k, p,  exp(1.0), alpha);

    // if(index == "PE" || index == "PC" || index == "XB")
    //   {ind_max = -ind;
    //   }else{
    //     ind_max = ind;
    //   }

    if(index == "PE" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }
  }else{

    ind = arma::datum::nan;
    warning("Some membership degrees are NaN (Suggestion: Increase the number of starting points RS)");

  }
  return List::create(Rcpp::Named("H") = H_opt,
                      Rcpp::Named("U") = U_opt,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k);

}


// [[Rcpp::export]]

List mainFKM_pf_noise_U(arma::mat data,
                        double b,
                        double delta,
                        unsigned int n,
                        unsigned int p,
                        unsigned int k,
                        arma::mat U,
                        double conv,
                        unsigned int maxit,
                        std::string index,
                        double alpha){

  int iter = 0;
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U_old = U;

  arma::rowvec d(n); d.zeros();
  unsigned int ki = 0;
  unsigned int kok = 0;
  arma::vec Uout(n); Uout.zeros();
  arma::uvec id;
  arma::uword mm = 0;

  arma::rowvec deltatemp(1); deltatemp = delta;

  double deltasq = pow(delta,2.0);
  double ind = 0;
  double ind_max = 0;

  bool prova = true;

  double value = 0;

  int nan_check = 1;

  while(prova && (iter < maxit))
  {

    iter++;
    U_old = U;

    H = centroids_FKM_pf(data, U, n, k, p, b);

    D = euclidean_distance(data, H, n, k);

    for(int i=0; i<n;i++)
    {

      nan_check = arma::is_finite(U.row(i));
      if(nan_check == 0){
        break;
      }
      d = D.row(i);

      if(d.min() == 0)
      {

        mm = d.index_min();
        U.row(i).zeros();
        U(i,mm) = 1;

      }else{

        ki = 0;
        kok = 0;
        d = arma::join_rows(d,deltatemp);
        d = sort(d);

        while(ki < k)

        {
          //if((d(ki) * (sum(1/d.subvec(0,ki)) + (1/deltasq))) <= ((1/b) + ki))
          if((d(ki) * (sum(1/d.subvec(0,ki)))) <= ((1/b) + ki))
          {
            kok = ki;
            ki++;
          }else{

            ki++;
          }
        }
        U.row(i) = 1/(1-b)*(((1+b*(kok))/(D.row(i)*(sum(1/d.subvec(0,kok)))))-b);
        U.row(i) = replace(U.row(i),0,0);

      }
      Uout(i) = 1 - sum(U.row(i));
    }

    prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

  }


  if(nan_check == 0)
  {
    ind = arma::datum::nan;
    warning("Some membership degrees are NaN (Suggestion: Increase the number of starting points RS)");


  }else{

    value = accu((((1-b)/(1+b))*pow(U,2.0) + ((2*b)/(1+b))*U) % D) + accu((((1-b)/(1+b))*pow(U,2.0)* deltasq)) + sum(((2*b)/(1+b))*Uout * deltasq);
    ind = indices(index, data, U, H, 1,  n, k, p,  exp(1.0), alpha);

    // if(index == "PE" || index == "PC" || index == "XB")
    //   {ind_max = -ind;
    //   }else{
    //     ind_max = ind;
    //   }

    if(index == "PE" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }
  }

  return List::create(Rcpp::Named("H") = H,
                      Rcpp::Named("U") = U,
                      Rcpp::Named("iter") = iter,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k);

}



///////////////////////

///////////////////////////////////
//          GKB FKM              //
///////////////////////////////////





// [[Rcpp::export]]

List mainFKM_gkb(arma::mat data,
                 double m,
                 double gam,
                 double mcn,
                 arma::vec vp,
                 unsigned int n,
                 unsigned int p,
                 unsigned int k,
                 unsigned rs,
                 double conv,
                 unsigned int maxit,
                 std::string index,
                 double alpha) {

  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U(n,k); U.zeros();
  arma::mat U_old(n,k); U_old.zeros();
  arma::cube F(p,p,k); F.zeros();

  double func = 0;
  double func_opt = 0;
  double ind = 0;
  double ind_max = 0;

  arma::mat H_opt(k, p); H.zeros();
  arma::mat U_opt(n,k); U.zeros();
  arma::cube F_opt(p,p,k); F_opt.zeros();

  arma::mat F0(p,p); F0.zeros();

  F0.diag().fill(pow(arma::det(arma::cov(data)), (1/(double)p)));

  for(int r=0; r<rs; r++)
  {



    int iter = 0;
    U = unifInit(n,k);
    U_old = U;

    bool prova = true;

    while(prova && (iter < maxit))
    {

      iter++;
      U_old = U;

      H = centroids_FKM(data, U, n, k, p, m);

      F = F_gkb(data, U,H, F0,  m,  gam, n, k, p, mcn, vp);

      D = euclidean_distance_gkb(data, H, F,n, k);

      U = memb_degree(D, m, n, k, p);

      prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

    }

    func = accu(pow(U,m)%D);
    it(r) = iter;
    value(r) = func;

    if(arma::is_finite(func) == true)
    {
      if ((r == 0) | (func < func_opt))
      {
        U_opt=U;
        H_opt=H;
        F_opt=F;
        func_opt=func;
      }
    }
  }
  ind = indices(index, data, U_opt, H_opt, m,  n, k, p,  exp(1.0), alpha);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }

  return List::create(Rcpp::Named("H") = H_opt,
                      Rcpp::Named("U") = U_opt,
                      Rcpp::Named("F") = F_opt,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("vp") = vp);

}

// [[Rcpp::export]]

List mainFKM_gkb_U(arma::mat data,
                   double m,
                   double gam,
                   double mcn,
                   arma::vec vp,
                   unsigned int n,
                   unsigned int p,
                   unsigned int k,
                   arma::mat U,
                   double conv,
                   unsigned int maxit,
                   std::string index,
                   double alpha) {

  int iter = 0;
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::cube F(p,p,k); F.zeros();
  arma::mat U_old = U;


  arma::mat F0(p,p); F0.zeros();

  F0.diag().fill(pow(arma::det(arma::cov(data)), (1/(double)p)));

  bool prova = true;

  double value = 0;
  double ind = 0;
  double ind_max = 0;

  while(prova && (iter < maxit))
  {

    iter++;
    U_old = U;

    H = centroids_FKM(data, U, n, k, p, m);

    F = F_gkb(data, U,H, F0,  m,  gam, n, k, p, mcn, vp);

    D = euclidean_distance_gkb(data, H, F,n, k);

    U = memb_degree(D, m, n, k, p);

    prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

  }

  value = accu(pow(U,m)%D);

  ind = indices(index, data, U, H, m,  n, k, p,  exp(1.0), alpha);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }
  return List::create(Rcpp::Named("H") = H,
                      Rcpp::Named("U") = U,
                      Rcpp::Named("F") = F,
                      Rcpp::Named("iter") = iter,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("vp") = vp);

}

///////////////////////

///////////////////////////////////
//  GKB ENTROPIC FKM              //
///////////////////////////////////




// [[Rcpp::export]]

List mainFKM_gkb_ent(arma::mat data,
                 double ent,
                 double gam,
                 double mcn,
                 arma::vec vp,
                 unsigned int n,
                 unsigned int p,
                 unsigned int k,
                 unsigned rs,
                 double conv,
                 unsigned int maxit,
                 std::string index,
                 double alpha) {

  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U(n,k); U.zeros();
  arma::mat U_old(n,k); U_old.zeros();
  arma::cube F(p,p,k); F.zeros();

  double func = 0;
  double func_opt = 0;
  double ind = 0;
  double ind_max = 0;

  arma::mat H_opt(k, p); H.zeros();
  arma::mat U_opt(n,k); U.zeros();
  arma::cube F_opt(p,p,k); F_opt.zeros();

  arma::mat F0(p,p); F0.zeros();

  F0.diag().fill(pow(arma::det(arma::cov(data)), (1/(double)p)));

  for(int r=0; r<rs; r++)
  {



    int iter = 0;
    U = unifInit(n,k);
    U_old = U;

    bool prova = true;

    while(prova && (iter < maxit))
    {

      iter++;
      U_old = U;

      H = centroids_FKM_ent(data, U, n, k, p);

      F = F_gkb_ent(data, U,H, F0,  gam, n, k, p, mcn, vp);

      D = euclidean_distance_gkb(data, H, F,n, k);

      U = memb_degree_ent(D, ent, n, k, p);

      prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

    }

    func = accu(U%D) + ent * accu(U % log(U));
    it(r) = iter;
    value(r) = func;

    if(arma::is_finite(func) == true)
    {
      if ((r == 0) | (func < func_opt))
      {
        U_opt=U;
        H_opt=H;
        F_opt=F;
        func_opt=func;
      }
    }
  }


  ind = indices(index, data, U_opt, H_opt, 2,  n, k, p,  exp(1.0), alpha);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }
  return List::create(Rcpp::Named("H") = H_opt,
                      Rcpp::Named("U") = U_opt,
                      Rcpp::Named("F") = F_opt,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("vp") = vp);

}

// [[Rcpp::export]]

List mainFKM_gkb_ent_U(arma::mat data,
                   double gam,
                   double mcn,
                   double ent,
                   arma::vec vp,
                   unsigned int n,
                   unsigned int p,
                   unsigned int k,
                   arma::mat U,
                   double conv,
                   unsigned int maxit,
                   std::string index,
                   double alpha) {

  int iter = 0;
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::cube F(p,p,k); F.zeros();
  arma::mat U_old = U;


  arma::mat F0(p,p); F0.zeros();

  F0.diag().fill(pow(arma::det(arma::cov(data)), (1/(double)p)));

  bool prova = true;

  double value = 0;
  double ind = 0;
  double ind_max = 0;

  while(prova && (iter < maxit))
  {

    iter++;
    U_old = U;

    H = centroids_FKM_ent(data, U, n, k, p);

    F = F_gkb_ent(data, U,H, F0,  gam, n, k, p, mcn, vp);

    D = euclidean_distance_gkb(data, H, F,n, k);

    U = memb_degree_ent(D, ent, n, k, p);

    prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

  }

  value = accu(U%D) + ent * accu(U % log(U));

  ind = indices(index, data, U, H, 2,  n, k, p,  exp(1.0), alpha);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }
  return List::create(Rcpp::Named("H") = H,
                      Rcpp::Named("U") = U,
                      Rcpp::Named("F") = F,
                      Rcpp::Named("iter") = iter,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("vp") = vp);

}


///////////////////////

///////////////////////////////////
//   GKB NOISE FKM              //
///////////////////////////////////

// [[Rcpp::export]]

List mainFKM_gkb_noise(arma::mat data,
                       double m,
                     double delta,
                     double gam,
                     double mcn,
                     arma::vec vp,
                     unsigned int n,
                     unsigned int p,
                     unsigned int k,
                     unsigned rs,
                     double conv,
                     unsigned int maxit,
                     std::string index,
                     double alpha) {

  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U(n,k); U.zeros();
  arma::mat U_old(n,k); U_old.zeros();
  arma::cube F(p,p,k); F.zeros();

  double func = 0;
  double func_opt = 0;
  double ind = 0;
  double ind_max = 0;

  arma::mat H_opt(k, p); H.zeros();
  arma::mat U_opt(n,k); U.zeros();
  arma::cube F_opt(p,p,k); F_opt.zeros();

  arma::mat F0(p,p); F0.zeros();

  F0.diag().fill(pow(arma::det(arma::cov(data)), (1/(double)p)));

  int nan_check = 1;

  for(int r=0; r<rs; r++)
  {
    int iter = 0;
    U = unifInit(n,k);
    U_old = U;

    arma::rowvec d(n); d.zeros();
    arma::vec Uout(n); Uout.zeros();
    arma::uword mm = 0;

    double q1 = 0;
    double q2 = 0;

    double deltasq = pow(delta,2.0);

    bool prova = true;

    while(prova && (iter < maxit))
    {

      iter++;
      U_old = U;

      H = centroids_FKM(data, U, n, k, p, m);

      F = F_gkb(data, U,H, F0, m, gam, n, k, p, mcn, vp);

      D = euclidean_distance_gkb(data, H, F,n, k);

      for(int i=0; i<n;i++)
      {

        d = D.row(i);

        if(d.min() == 0)
        {

          mm = d.index_min();
          U.row(i).zeros();
          U(i,mm) = 1;

        }else{

          for(int j=0; j<k;j++)
          {

            q1 = pow(1/arma::as_scalar(D(i,j)),1/(m-1.0)) / sum(pow(1/D.row(i),1/(m-1.0)));
            q2 = pow(arma::as_scalar(D(i,j)),1/(m-1.0)) / pow(delta,2/(m-1.0));
            U(i,j) = 1/((1/q1)+q2);

          }
        }

        Uout(i) = 1 - sum(U.row(i));
        Uout(i) = pow(Uout(i),m) * deltasq;

      }

      prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

    }

    if(nan_check == 0)
    {
      func = arma::datum::nan;
    }else{

      func = accu(pow(U,m)%D) + sum(Uout);
    }

    it(r) = iter;
    value(r) = func;

    if(arma::is_finite(func) == true)
    {
      if ((r == 0) | (func < func_opt))
      {
        U_opt=U;
        H_opt=H;
        F_opt=F;
        func_opt=func;
      }
    }
  }

  if(arma::is_finite(func_opt)){

    ind = indices(index, data, U_opt, H_opt, m,  n, k, p,  exp(1.0), alpha);

    // if(index == "PE" || index == "PC" || index == "XB")
    //   {ind_max = -ind;
    //   }else{
    //     ind_max = ind;
    //   }

    if(index == "PE" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }
  }else{

    ind = arma::datum::nan;
    warning("Some membership degrees are NaN (Suggestion: Increase the number of starting points RS)");

  }


  return List::create(Rcpp::Named("H") = H_opt,
                      Rcpp::Named("U") = U_opt,
                      Rcpp::Named("F") = F_opt,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("vp") = vp);


}

// [[Rcpp::export]]

List mainFKM_gkb_noise_U(arma::mat data,
                         double m,
                       double gam,
                       double mcn,
                       double delta,
                       arma::vec vp,
                       unsigned int n,
                       unsigned int p,
                       unsigned int k,
                       arma::mat U,
                       double conv,
                       unsigned int maxit,
                       std::string index,
                       double alpha) {

  int iter = 0;
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::cube F(p,p,k); F.zeros();
  arma::mat U_old = U;


  arma::mat F0(p,p); F0.zeros();

  F0.diag().fill(pow(arma::det(arma::cov(data)), (1/(double)p)));

  arma::rowvec d(n); d.zeros();
  arma::vec Uout(n); Uout.zeros();
  double q1 = 0;
  double q2 = 0;

  double deltasq = pow(delta,2.0);

  bool prova = true;

  double value = 0;
  double ind = 0;
  double ind_max = 0;

  arma::uword mm = 0;

  int nan_check = 1;

  while(prova && (iter < maxit))
  {

    iter++;
    U_old = U;

    H = centroids_FKM(data, U, n, k, p, m);

    F = F_gkb(data, U,H, F0, m,  gam, n, k, p, mcn, vp);

    D = euclidean_distance_gkb(data, H, F,n, k);

    for(int i=0; i<n;i++)
    {

      d = D.row(i);

      if(d.min() == 0)
      {

        mm = d.index_min();
        U.row(i).zeros();
        U(i,mm) = 1;

      }else{

        for(int j=0; j<k;j++)
        {

          q1 = pow(1/arma::as_scalar(D(i,j)),1/(m-1.0)) / sum(pow(1/D.row(i),1/(m-1.0)));
          q2 = pow(arma::as_scalar(D(i,j)),1/(m-1.0)) / pow(delta,2/(m-1.0));
          nan_check = arma::is_finite(U.row(i));
          if(nan_check == 0){
            break;
          }
          U(i,j) = 1/((1/q1)+q2);

        }
      }

      Uout(i) = 1 - sum(U.row(i));
      Uout(i) = pow(Uout(i),m) * deltasq;

    }

    prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

  }

  if(nan_check == 0)
  {
    ind = arma::datum::nan;
    warning("Some membership degrees are NaN (Suggestion: Increase the number of starting points RS)");


  }else{

    value = accu(pow(U,m)%D) + sum(Uout);
    ind = indices(index, data, U, H, m,  n, k, p,  exp(1.0), alpha);

    if(index == "PE" || index == "PC" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }
  }

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }
  return List::create(Rcpp::Named("H") = H,
                      Rcpp::Named("U") = U,
                      Rcpp::Named("F") = F,
                      Rcpp::Named("iter") = iter,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("vp") = vp);



}

///////////////////////

///////////////////////////////////
//GKB NOISE-ENTROPIC FKM          //
///////////////////////////////////

// [[Rcpp::export]]

List mainFKM_gkb_ent_noise(arma::mat data,
                       double ent,
                       double delta,
                       double gam,
                       double mcn,
                       arma::vec vp,
                       unsigned int n,
                       unsigned int p,
                       unsigned int k,
                       unsigned rs,
                       double conv,
                       unsigned int maxit,
                       std::string index,
                       double alpha) {

  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U(n,k); U.zeros();
  arma::mat U_old(n,k); U_old.zeros();
  arma::cube F(p,p,k); F.zeros();

  double func = 0;
  double func_opt = 0;
  double ind = 0;
  double ind_max = 0;

  arma::mat H_opt(k, p); H.zeros();
  arma::mat U_opt(n,k); U.zeros();
  arma::cube F_opt(p,p,k); F_opt.zeros();

  arma::mat F0(p,p); F0.zeros();

  arma::rowvec d(n); d.zeros();
  arma::vec Uout(n); Uout.zeros();
  arma::uword mm = 0;

  F0.diag().fill(pow(arma::det(arma::cov(data)), (1/(double)p)));

  for(int r=0; r<rs; r++)
  {
    int iter = 0;
    U = unifInit(n,k);
    U_old = U;

    double eps = std::numeric_limits<double>::epsilon();

    bool prova = true;

    while(prova && (iter < maxit))
    {

      iter++;
      U_old = U;

      H = centroids_FKM_ent(data, U, n, k, p);

      F = F_gkb_ent(data, U,H, F0, gam, n, k, p, mcn, vp);

      D = euclidean_distance_gkb(data, H, F,n, k);

      for(int i=0; i<n;i++)
      {

        d = D.row(i);

        if(d.min() == 0)
        {

          mm = d.index_min();
          U.row(i).zeros();
          U(i,mm) = 1;

        }else{

          for(int j=0; j<k;j++)
          {
            U(i,j) = exp(-arma::as_scalar(D(i,j))/ent) / (sum(exp(D.row(i)/(-ent))) + exp(- pow(delta,2.0) / ent));

            if(arma::is_finite(U(i,j)) == false)
            {
              stop("Some membership degrees are NaN (Suggestion: run FKM.gkb.ent.noise using standardized data)");
            }

            if(U(i,j) < eps)
            {
              U(i,j) = eps;
            }
          }
          Uout(i) = 1 - sum(U.row(i));
        }
      }

      prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

    }

    func = accu(U%D) + ent * accu(U % log(U)) + sum(Uout) * pow(delta,2.0);
    it(r) = iter;
    value(r) = func;

    if(arma::is_finite(func) == true)
    {
      if ((r == 0) | (func < func_opt))
      {
        U_opt=U;
        H_opt=H;
        F_opt=F;
        func_opt=func;
      }
    }



  }


  ind = indices(index, data, U_opt, H_opt, 2,  n, k, p,  exp(1.0), alpha);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }
  return List::create(Rcpp::Named("H") = H_opt,
                      Rcpp::Named("U") = U_opt,
                      Rcpp::Named("F") = F_opt,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("vp") = vp);

}

// [[Rcpp::export]]

List mainFKM_gkb_ent_noise_U(arma::mat data,
                         double ent,
                         double delta,
                         double gam,
                         double mcn,
                         arma::vec vp,
                         unsigned int n,
                         unsigned int p,
                         unsigned int k,
                         arma::mat U,
                         double conv,
                         unsigned int maxit,
                         std::string index,
                         double alpha) {

  int iter = 0;
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::cube F(p,p,k); F.zeros();
  arma::mat U_old = U;
  arma::uword mm = 0;


  arma::mat F0(p,p); F0.zeros();

  F0.diag().fill(pow(arma::det(arma::cov(data)), (1/(double)p)));

  arma::rowvec d(n); d.zeros();
  arma::vec Uout(n); Uout.zeros();

  bool prova = true;

  double eps = std::numeric_limits<double>::epsilon();

  double value = 0;
  double ind = 0;
  double ind_max = 0;

  while(prova && (iter < maxit))
  {

    iter++;
    U_old = U;

    H = centroids_FKM_ent(data, U, n, k, p);

    F = F_gkb_ent(data, U,H, F0, gam, n, k, p, mcn, vp);

    D = euclidean_distance_gkb(data, H, F,n, k);

    for(int i=0; i<n;i++)
    {

      d = D.row(i);

      if(d.min() == 0)
      {

        mm = d.index_min();
        U.row(i).zeros();
        U(i,mm) = 1;

      }else{

        for(int j=0; j<k;j++)
        {
          U(i,j) = exp(-arma::as_scalar(D(i,j))/ent) / (sum(exp(D.row(i)/(-ent))) + exp(- pow(delta,2.0) / ent));

          if(arma::is_finite(U(i,j)) == false)
          {
            stop("Some membership degrees are NaN (Suggestion: run FKM.gkb.ent.noise using standardized data)");
          }

          if(U(i,j) < eps)
          {
            U(i,j) = eps;
          }
        }
        Uout(i) = 1 - sum(U.row(i));
      }
    }

    prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

  }

  value = accu(U%D) + ent * accu(U % log(U)) + sum(Uout) * pow(delta,2.0);

  ind = indices(index, data, U, H, 2,  n, k, p,  exp(1.0), alpha);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }
  return List::create(Rcpp::Named("H") = H,
                      Rcpp::Named("U") = U,
                      Rcpp::Named("F") = F,
                      Rcpp::Named("iter") = iter,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("vp") = vp);



}



/////////////////////////////////////////////////


///////////////////////////////////
//          FKM GK              //
///////////////////////////////////




// [[Rcpp::export]]

List mainFKM_gk(arma::mat data,
                double m,
                arma::vec vp,
                unsigned int n,
                unsigned int p,
                unsigned int k,
                unsigned rs,
                double conv,
                unsigned int maxit,
                std::string index,
                double alpha) {

  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat D_old(n, k); D_old.zeros();
  arma::mat U(n,k); U.zeros();
  arma::mat U_old(n,k); U_old.zeros();
  arma::cube F(p,p,k); F.zeros();

  arma::mat D_temp(n,k); D_temp.zeros();

  double func = 0;
  double func_opt = 0;
  double ind = 0;
  double ind_max = 0;

  arma::mat H_opt(k, p); H.zeros();
  arma::mat U_opt(n,k); U.zeros();
  arma::cube F_opt(p,p,k); F_opt.zeros();

  arma::cube F_old(p,p,k); F_old.zeros();
  arma::mat H_old(k, p); H_old.zeros();

  int countNoConv = 0;
  int warn = 0;

  for(int r=0; r<rs; r++)
  {



    int iter = 0;
    U = unifInit(n,k);
    U_old = U;

    bool prova = true;

    while(prova && (iter < maxit))
    {

      iter++;
      U_old = U;
      F_old = F;
      H_old = H;
      D_old = D;

      H = centroids_FKM(data, U, n, k, p, m);

      F = F_gk(data,U, H, m,n, k, p, vp);

      D_temp = euclidean_distance_gk(data, H, F,D, n, k,p);

      if(D_temp.is_empty()){

        U = U_old;
        D = D_old;
        F = F_old;
        H = H_old;
        iter = iter - 1;
        countNoConv++;

      }else{

        D = D_temp;
        U = memb_degree(D, m, n, k, p);

      }

      prova = arma::as_scalar(accu(abs(U_old - U)))>conv;


    }

    if(iter == 1)
    {
      func = arma::datum::nan;
    }else{
      func = accu(pow(U,m)%D);
    }

    it(r) = iter;
    value(r) = func;

    if(arma::is_finite(func) == true)
    {
      if ((r == 0) | (func < func_opt))
      {
        U_opt=U;
        H_opt=H;
        F_opt=F;
        func_opt=func;
      }
    }
  }

  if(arma::is_finite(func_opt)){

    ind = indices(index, data, U_opt, H_opt, m,  n, k, p,  exp(1.0), alpha);

    // if(index == "PE" || index == "PC" || index == "XB")
    //   {ind_max = -ind;
    //   }else{
    //     ind_max = ind;
    //   }

    if(index == "PE" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }
  }else{

    ind = arma::datum::nan;
    warning("Some membership degrees are NaN (Suggestion: Increase the number of starting points RS)");
  }

  if(countNoConv == rs)
  {
    warn = 1;
  }

  return List::create(Rcpp::Named("H") = H_opt,
                      Rcpp::Named("U") = U_opt,
                      Rcpp::Named("F") = F_opt,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("vp") = vp,
                      Rcpp::Named("warn") = warn);

}

// [[Rcpp::export]]

List mainFKM_gk_U(arma::mat data,
                  double m,
                  arma::vec vp,
                  unsigned int n,
                  unsigned int p,
                  unsigned int k,
                  arma::mat U,
                  double conv,
                  unsigned int maxit,
                  std::string index,
                  double alpha) {

  int iter = 0;
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat D_old(n, k); D_old.zeros();
  arma::cube F(p,p,k); F.zeros();
  arma::mat U_old = U;

  arma::mat D_temp(n,k); D_temp.zeros();
  arma::cube F_old(p,p,k); F_old.zeros();
  arma::mat H_old(k, p); H_old.zeros();

  bool prova = true;

  double value = 0;
  double ind = 0;
  double ind_max = 0;

  int warn = 0;

  while(prova && (iter < maxit))
  {

    iter++;
    U_old = U;
    F_old = F;
    H_old = H;
    D_old = D;

    H = centroids_FKM(data, U, n, k, p, m);

    F = F_gk(data,U, H,m, n, k, p, vp);

    D_temp = euclidean_distance_gk(data, H, F,D, n, k,p);

    if(D_temp.is_empty()){

      U = U_old;
      D = D_old;
      warn = 1;

    }else{

      D = D_temp;
      U = memb_degree(D, m, n, k, p);

    }

    prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

  }

  if(iter == 1)
  {
    value = arma::datum::nan;
    ind = arma::datum::nan;
  }else{
    value = accu(pow(U,m)%D);
    ind = indices(index, data, U, H, m,  n, k, p,  exp(1.0), alpha);

    // if(index == "PE" || index == "PC" || index == "XB")
    //   {ind_max = -ind;
    //   }else{
    //     ind_max = ind;
    //   }

    if(index == "PE" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }
  }


  return List::create(Rcpp::Named("H") = H,
                      Rcpp::Named("U") = U,
                      Rcpp::Named("F") = F,
                      Rcpp::Named("iter") = iter,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("vp") = vp,
                      Rcpp::Named("warn") = warn);

}


///////////////////////


///////////////////////////////////
//    FKM GK ENTROPIC            //
///////////////////////////////////




// [[Rcpp::export]]

List mainFKM_gk_ent(arma::mat data,
                double ent,
                arma::vec vp,
                unsigned int n,
                unsigned int p,
                unsigned int k,
                unsigned rs,
                double conv,
                unsigned int maxit,
                std::string index,
                double alpha) {

  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat D_old(n, k); D_old.zeros();
  arma::mat U(n,k); U.zeros();
  arma::mat U_old(n,k); U_old.zeros();
  arma::cube F(p,p,k); F.zeros();

  arma::mat D_temp(n,k); D_temp.zeros();

  double func = 0;
  double func_opt = 0;
  double ind = 0;
  double ind_max = 0;

  arma::mat H_opt(k, p); H.zeros();
  arma::mat U_opt(n,k); U.zeros();
  arma::cube F_opt(p,p,k); F_opt.zeros();

  arma::cube F_old(p,p,k); F_old.zeros();
  arma::mat H_old(k, p); H_old.zeros();

  int countNoConv = 0;
  int warn = 0;

  for(int r=0; r<rs; r++)
  {

    int iter = 0;
    U = unifInit(n,k);
    U_old = U;

    bool prova = true;

    while(prova && (iter < maxit))
    {

      iter++;
      U_old = U;
      D_old = D;


      H = centroids_FKM_ent(data, U, n, k, p);

      F = F_gk_ent(data,U, H,n, k, p, vp);

      D_temp = euclidean_distance_gk(data, H, F,D, n, k,p);

      if(D_temp.is_empty()){

        U = U_old;
        D = D_old;
        F = F_old;
        H = H_old;
        iter = iter - 1;
        countNoConv++;

      }else{

        D = D_temp;
        U = memb_degree_ent(D, ent,n, k, p);

      }
      prova = arma::as_scalar(accu(abs(U_old - U)))>conv;
    }

    if(iter == 1)
    {
      func = arma::datum::nan;
    }else{
      func =  accu(U%D) + ent * accu(U % log(U));
    }

    it(r) = iter;
    value(r) = func;

    if(arma::is_finite(func) == true)
    {
      if ((r == 0) | (func < func_opt))
      {
        U_opt=U;
        H_opt=H;
        F_opt=F;
        func_opt=func;
      }
    }
  }


  if(arma::is_finite(func_opt)){

    ind = indices(index, data, U_opt, H_opt, 1,  n, k, p,  exp(1.0), alpha);

    // if(index == "PE" || index == "PC" || index == "XB")
    //   {ind_max = -ind;
    //   }else{
    //     ind_max = ind;
    //   }

    if(index == "PE" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }
  }else{

    ind = arma::datum::nan;
    warning("Some membership degrees are NaN (Suggestion: Increase the number of starting points RS)");

  }

  if(countNoConv == rs)
  {
    warn = 1;
  }

  return List::create(Rcpp::Named("H") = H_opt,
                      Rcpp::Named("U") = U_opt,
                      Rcpp::Named("F") = F_opt,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("vp") = vp,
                      Rcpp::Named("warn") = warn);

}

// [[Rcpp::export]]

List mainFKM_gk_ent_U(arma::mat data,
                  double ent,
                  arma::vec vp,
                  unsigned int n,
                  unsigned int p,
                  unsigned int k,
                  arma::mat U,
                  double conv,
                  unsigned int maxit,
                  std::string index,
                  double alpha) {

  int iter = 0;
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat D_old(n, k); D_old.zeros();
  arma::cube F(p,p,k); F.zeros();
  arma::mat U_old = U;

  arma::mat D_temp(n,k); D_temp.zeros();

  bool prova = true;

  double value = 0;
  double ind = 0;
  double ind_max = 0;

  int warn = 0;

  while(prova && (iter < maxit))
  {

    iter++;
    U_old = U;
    D_old = D;

    H = centroids_FKM_ent(data, U, n, k, p);

    F = F_gk_ent(data,U, H,n, k, p, vp);

    D_temp = euclidean_distance_gk(data, H, F,D, n, k,p);

    if(D_temp.is_empty()){

      U = U_old;
      D = D_old;
      warn = 1;

    }else{

      D = D_temp;
      U = memb_degree_ent(D, ent,n, k, p);
    }
    prova = arma::as_scalar(accu(abs(U_old - U)))>conv;
  }

  if(iter == 1)
  {
    value = arma::datum::nan;
    ind = arma::datum::nan;
  }else{
    value = accu(U%D) + ent * accu(U % log(U));
    ind = indices(index, data, U, H, 2,  n, k, p,  exp(1.0), alpha);

    // if(index == "PE" || index == "PC" || index == "XB")
    //   {ind_max = -ind;
    //   }else{
    //     ind_max = ind;
    //   }

    if(index == "PE" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }
  }

  return List::create(Rcpp::Named("H") = H,
                      Rcpp::Named("U") = U,
                      Rcpp::Named("F") = F,
                      Rcpp::Named("iter") = iter,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("vp") = vp,
                      Rcpp::Named("warn") = warn);
}


///////////////////////



///////////////////////////////////
//      FKM GK NOISE            //
///////////////////////////////////

// [[Rcpp::export]]

List mainFKM_gk_noise(arma::mat data,
                    double m,
                    double delta,
                    arma::vec vp,
                    unsigned int n,
                    unsigned int p,
                    unsigned int k,
                    unsigned rs,
                    double conv,
                    unsigned int maxit,
                    std::string index,
                    double alpha) {

  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat D_old(n, k); D_old.zeros();
  arma::mat U(n,k); U.zeros();
  arma::cube F(p,p,k); F.zeros();

  arma::mat D_temp(n,k); D_temp.zeros();

  arma::rowvec d(n); d.zeros();
  arma::vec Uout(n); Uout.zeros();
  double q1 = 0;
  double q2 = 0;

  double deltasq = pow(delta,2.0);

  double func = 0;
  double func_opt = 0;
  double ind = 0;
  double ind_max = 0;

  arma::mat H_opt(k, p); H.zeros();
  arma::mat U_opt(n,k); U.zeros();
  arma::mat U_old(n,k); U_old.zeros();
  arma::cube F_opt(p,p,k); F_opt.zeros();

  arma::cube F_old(p,p,k); F_old.zeros();
  arma::mat H_old(k, p); H_old.zeros();

  int countNoConv = 0;
  int warn = 0;

  int nan_check = 1;


  for(int r=0; r<rs; r++)
  {

    int iter = 0;
    U = unifInit(n,k);
    U_old = U;

    bool prova = true;

    while(prova && (iter < maxit))
    {

      iter++;
      U_old = U;
      D_old = D;

      H = centroids_FKM(data, U, n, k, p, m);

      F = F_gk(data,U, H, m,n, k, p, vp);

      D_temp = euclidean_distance_gk(data, H, F,D, n, k,p);

      if(D_temp.is_empty()){

        U = U_old;
        D = D_old;
        F = F_old;
        H = H_old;
        iter = iter - 1;
        countNoConv++;

      }else{

        D = D_temp;
        for(int i=0; i<n;i++)
        {

          d = D.row(i);

          if(d.min() == 0)
          {

            arma::uword mm = d.index_min();
            U.row(i).zeros();
            U(i,mm) = 1;

          }else{

            for(int j=0; j<k;j++)
            {

              q1 = pow(1/arma::as_scalar(D(i,j)),1/(m-1.0)) / sum(pow(1/D.row(i),1/(m-1.0)));
              q2 = pow(arma::as_scalar(D(i,j)),1/(m-1.0)) / pow(delta,2/(m-1.0));
              U(i,j) = 1/((1/q1)+q2);

            }
          }

          Uout(i) = 1 - sum(U.row(i));
          Uout(i) = pow(Uout(i),m) * deltasq;
        }
      }
      prova = arma::as_scalar(accu(abs(U_old - U)))>conv;
    }

    if((iter == 1) | (nan_check == 0))
    {
      func = arma::datum::nan;
    }else{
      func =  accu(pow(U,m)%D) + sum(Uout);
    }

    it(r) = iter;
    value(r) = func;

    if(arma::is_finite(func) == true)
    {
      if ((r == 0) | (func < func_opt))
      {
        U_opt=U;
        H_opt=H;
        F_opt=F;
        func_opt=func;
      }
    }
  }


  if(arma::is_finite(func_opt)){

    ind = indices(index, data, U_opt, H_opt, m,  n, k, p,  exp(1.0), alpha);

    // if(index == "PE" || index == "PC" || index == "XB")
    //   {ind_max = -ind;
    //   }else{
    //     ind_max = ind;
    //   }

    if(index == "PE" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }
  }else{

    ind = arma::datum::nan;
    warning("Some membership degrees are NaN (Suggestion: Increase the number of starting points RS)");

  }
  if(countNoConv == rs)
  {
    warn = 1;
  }
  return List::create(Rcpp::Named("H") = H_opt,
                      Rcpp::Named("U") = U_opt,
                      Rcpp::Named("F") = F_opt,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("vp") = vp,
                      Rcpp::Named("warn") = warn);

}

// [[Rcpp::export]]

List mainFKM_gk_noise_U(arma::mat data,
                      double m,
                      double delta,
                      arma::vec vp,
                      unsigned int n,
                      unsigned int p,
                      unsigned int k,
                      arma::mat U,
                      double conv,
                      unsigned int maxit,
                      std::string index,
                      double alpha) {

  int iter = 0;
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat D_old(n, k); D_old.zeros();
  arma::cube F(p,p,k); F.zeros();
  arma::mat U_old = U;

  arma::mat D_temp(n,k); D_temp.zeros();

  bool prova = true;

  double value = 0;
  double ind = 0;
  double ind_max = 0;

  arma::rowvec d(n); d.zeros();
  arma::vec Uout(n); Uout.zeros();
  double q1 = 0;
  double q2 = 0;

  double deltasq = pow(delta,2.0);

  int nan_check = 1;

  int warn = 0;

  while(prova && (iter < maxit))
  {

    iter++;
    U_old = U;
    D_old = D;

    H = centroids_FKM(data, U, n, k, p, m);

    F = F_gk(data,U, H, m,n, k, p, vp);

    D_temp = euclidean_distance_gk(data, H, F,D, n, k,p);

    if(D_temp.is_empty()){

      U = U_old;
      D = D_old;
      warn = 1;

    }else{

      D = D_temp;
      for(int i=0; i<n;i++)
      {

        d = D.row(i);
        if(d.min() == 0)
        {

          arma::uword mm = d.index_min();
          U.row(i).zeros();
          U(i,mm) = 1;

        }else{

          for(int j=0; j<k;j++)
          {

            q1 = pow(1/arma::as_scalar(D(i,j)),1/(m-1.0)) / sum(pow(1/D.row(i),1/(m-1.0)));
            q2 = pow(arma::as_scalar(D(i,j)),1/(m-1.0)) / pow(delta,2/(m-1.0));

            nan_check = arma::is_finite(U.row(i));
            if(nan_check == 0){
              break;
            }

            U(i,j) = 1/((1/q1)+q2);

          }
        }

        Uout(i) = 1 - sum(U.row(i));
        Uout(i) = pow(Uout(i),m) * deltasq;
      }
    }
    prova = arma::as_scalar(accu(abs(U_old - U)))>conv;
  }

  if( (iter == 1) | (nan_check == 0))
  {
    value = arma::datum::nan;
    ind = arma::datum::nan;
  }else{
    value = accu(pow(U,m)%D) + sum(Uout);
    ind = indices(index, data, U, H, m,  n, k, p,  exp(1.0), alpha);

    // if(index == "PE" || index == "PC" || index == "XB")
    //   {ind_max = -ind;
    //   }else{
    //     ind_max = ind;
    //   }

    if(index == "PE" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }
  }


  return List::create(Rcpp::Named("H") = H,
                      Rcpp::Named("U") = U,
                      Rcpp::Named("F") = F,
                      Rcpp::Named("iter") = iter,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("vp") = vp,
                      Rcpp::Named("warn") = warn);

}




///////////////////////



///////////////////////////////////
//FKM GK NOISE - ENTROPIC        //
///////////////////////////////////


// [[Rcpp::export]]

List mainFKM_gk_ent_noise(arma::mat data,
                      double ent,
                      double delta,
                      arma::vec vp,
                      unsigned int n,
                      unsigned int p,
                      unsigned int k,
                      unsigned rs,
                      double conv,
                      unsigned int maxit,
                      std::string index,
                      double alpha) {

  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat D_old(n, k); D_old.zeros();
  arma::mat U(n,k); U.zeros();
  arma::cube F(p,p,k); F.zeros();

  arma::mat D_temp(n,k); D_temp.zeros();

  arma::rowvec d(n); d.zeros();
  arma::vec Uout(n); Uout.zeros();

  double func = 0;
  double func_opt = 0;
  double ind = 0;
  double ind_max = 0;

  arma::mat H_opt(k, p); H.zeros();
  arma::mat U_opt(n,k); U.zeros();
  arma::mat U_old(n,k); U_old.zeros();
  arma::cube F_opt(p,p,k); F_opt.zeros();

  arma::cube F_old(p,p,k); F_old.zeros();
  arma::mat H_old(k, p); H_old.zeros();

  arma::uword mm = 0;

  int countNoConv = 0;
  int warn = 0;

  double eps = std::numeric_limits<double>::epsilon();

  for(int r=0; r<rs; r++)
  {

    int iter = 0;
    U = unifInit(n,k);
    U_old = U;

    bool prova = true;

    while(prova && (iter < maxit))
    {

      iter++;
      U_old = U;
      D_old = D;


      H = centroids_FKM_ent(data, U, n, k, p);

      F = F_gk_ent(data,U, H,n, k, p, vp);

      D_temp = euclidean_distance_gk(data, H, F,D, n, k,p);

      if(D_temp.is_empty()){

        U = U_old;
        D = D_old;
        F = F_old;
        H = H_old;
        iter = iter - 1;
        countNoConv++;

      }else{

        D = D_temp;
        for(int i=0; i<n;i++)
        {

          d = D.row(i);

          if(d.min() == 0)
          {

            mm = d.index_min();
            U.row(i).zeros();
            U(i,mm) = 1;

          }else{

            for(int j=0; j<k;j++)
            {
              U(i,j) = exp(-arma::as_scalar(D(i,j))/ent) / (sum(exp(D.row(i)/(-ent))) + exp(- pow(delta,2.0) / ent));

              if(arma::is_finite(U(i,j)) == false)
              {
                stop("Some membership degrees are NaN (Suggestion: run FKM.gk.ent using standardized data)");
              }

              if(U(i,j) < eps)
              {
                U(i,j) = eps;
              }
            }
            Uout(i) = 1 - sum(U.row(i));
          }
        }
      }
      prova = arma::as_scalar(accu(abs(U_old - U)))>conv;
    }

    if(iter == 1 && D_temp.is_empty())
    {
      func = arma::datum::nan;
    }else{
      func =  accu(U%D) + ent * accu(U % log(U)) + sum(Uout) * pow(delta,2.0);
    }

    it(r) = iter;
    value(r) = func;

    if(arma::is_finite(func) == true)
    {
      if ((r == 0) | (func < func_opt))
      {
        U_opt=U;
        H_opt=H;
        F_opt=F;
        func_opt=func;
      }
    }
  }


  ind = indices(index, data, U_opt, H_opt, 2,  n, k, p,  exp(1.0), alpha);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }
  return List::create(Rcpp::Named("H") = H_opt,
                      Rcpp::Named("U") = U_opt,
                      Rcpp::Named("F") = F_opt,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("vp") = vp,
                      Rcpp::Named("warn") = warn);

}

// [[Rcpp::export]]

List mainFKM_gk_ent_noise_U(arma::mat data,
                        double ent,
                        double delta,
                        arma::vec vp,
                        unsigned int n,
                        unsigned int p,
                        unsigned int k,
                        arma::mat U,
                        double conv,
                        unsigned int maxit,
                        std::string index,
                        double alpha) {

  int iter = 0;
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat D_old(n, k); D_old.zeros();
  arma::cube F(p,p,k); F.zeros();
  arma::mat U_old = U;

  arma::mat D_temp(n,k); D_temp.zeros();

  bool prova = true;

  double value = 0;
  double ind = 0;
  double ind_max = 0;

  arma::rowvec d(n); d.zeros();
  arma::vec Uout(n); Uout.zeros();

  arma::uword mm = 0;

  double eps = std::numeric_limits<double>::epsilon();

  int warn = 0;

  while(prova && (iter < maxit))
  {

    iter++;
    U_old = U;
    D_old = D;

    H = centroids_FKM_ent(data, U, n, k, p);

    F = F_gk_ent(data,U, H,n, k, p, vp);

    D_temp = euclidean_distance_gk(data, H, F,D, n, k,p);

    if(D_temp.is_empty()){

      U = U_old;
      D = D_old;
      warn = 1;

    }else{

      D = D_temp;
      for(int i=0; i<n;i++)
      {

        d = D.row(i);

        if(d.min() == 0)
        {

          mm = d.index_min();
          U.row(i).zeros();
          U(i,mm) = 1;

        }else{

          for(int j=0; j<k;j++)
          {
            U(i,j) = exp(-arma::as_scalar(D(i,j))/ent) / (sum(exp(D.row(i)/(-ent))) + exp(- pow(delta,2.0) / ent));

            if(arma::is_finite(U(i,j)) == false)
            {
              stop("Some membership degrees are NaN (Suggestion: run FKM.gk.ent.noise using standardized data)");
            }

            if(U(i,j) < eps)
            {
              U(i,j) = eps;
            }
          }
          Uout(i) = 1 - sum(U.row(i));
        }
      }
    }
    prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

  }

  if(iter == 1 && D_temp.is_empty())
  {
    value = arma::datum::nan;
    ind = arma::datum::nan;
  }else{
    value = accu(U%D) + ent * accu(U % log(U)) + sum(Uout) * pow(delta,2.0);
    ind = indices(index, data, U, H, 2,  n, k, p,  exp(1.0), alpha);

    // if(index == "PE" || index == "PC" || index == "XB")
    //   {ind_max = -ind;
    //   }else{
    //     ind_max = ind;
    //   }

    if(index == "PE" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }
  }

  return List::create(Rcpp::Named("H") = H,
                      Rcpp::Named("U") = U,
                      Rcpp::Named("F") = F,
                      Rcpp::Named("iter") = iter,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("vp") = vp);

}



/////////////////////////////////////////////////



///////////////////////////////////
//         FUZZY  MEDOIDS        //
///////////////////////////////////


// [[Rcpp::export]]

List mainFKM_med(arma::mat data,
                 double m,
                 unsigned int n,
                 unsigned int p,
                 unsigned int k,
                 unsigned rs,
                 double conv,
                 unsigned int maxit,
                 std::string index,
                 double alpha) {

  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U(n,k); U.zeros();

  double func = 0;
  double func_opt = 0;

  arma::mat H_opt(k, p); H.zeros();
  arma::mat U_opt(n,k); U.zeros();
  arma::mat U_old(n,k); U_old.zeros();
  arma::uvec medoid_opt(k); medoid_opt.zeros();
  double min_med_const = pow(10.0,5.0) * accu(pow(data,2.0));
  double min_med_old = 0;
  double min_med = 0;

  arma::uvec medoid(k); medoid.zeros();

  bool index_temp_med = true;
  bool value_temp_med = true;

  double ind = 0;
  double ind_max = 0;

  int nan_check = 1;

  for(int r=0; r<rs; r++)
  {
    int iter = 0;
    U = unifInit(n,k);
    U_old = U;

    bool prova = true;

    while(prova && (iter < maxit))
    {

      iter++;
      U_old = U;

      for(int c = 0; c < k; c++)
      {
        min_med_old = min_med_const;

        for(int i=0; i<n; i++)
        {
          min_med = 0;
          nan_check = arma::is_finite(U.row(i));
          if(nan_check == 0){
            break;
          }

          for(int j=0; j<n; j++)
          {
            min_med = min_med + pow(U(j,c),m) * sum(pow(data.row(i) - data.row(j),2.0));
          }

          index_temp_med = Match(i,medoid);
          value_temp_med = (min_med < min_med_old);

          if(index_temp_med & value_temp_med)
          {
            min_med_old = min_med;
            medoid(c) = i;
          }
        }

        H.row(c) = data.row(medoid(c));
      }

      D = euclidean_distance_medoid(data, H, n, k);

      U = memb_degree_medoid(D, medoid, m, n, k, p);

      prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

    }

    if(nan_check == 0)
    {
      func = arma::datum::nan;
    }else{

      func = accu(pow(U,m)%D);
    }

    it(r) = iter;
    value(r) = func;

    if(arma::is_finite(func) == true)
    {
      if ((r == 0) | (func < func_opt))
      {
        U_opt=U;
        H_opt=H;
        medoid_opt = medoid;
        func_opt=func;
      }
    }
  }

  if(arma::is_finite(func_opt)){

    ind = indices(index, data, U_opt, H_opt, m,  n, k, p,  exp(1.0), alpha);

    // if(index == "PE" || index == "PC" || index == "XB")
    //   {ind_max = -ind;
    //   }else{
    //     ind_max = ind;
    //   }

    if(index == "PE" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }
  }else{

    ind = arma::datum::nan;
    warning("Some membership degrees are NaN (Suggestion: Increase the number of starting points RS)");

  }

  return List::create(Rcpp::Named("H") = H_opt,
                      Rcpp::Named("U") = U_opt,
                      Rcpp::Named("medoid") = medoid_opt + 1,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k);

}


// [[Rcpp::export]]

List mainFKM_med_U(arma::mat data,
                   double m,
                   unsigned int n,
                   unsigned int p,
                   unsigned int k,
                   arma::mat U,
                   double conv,
                   unsigned int maxit,
                   std::string index,
                   double alpha) {

  int iter = 0;
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U_old = U;

  bool prova = true;

  double value = 0;

  double min_med_const = pow(10.0,5.0) * accu(pow(data,2.0));
  double min_med_old = 0;
  double min_med = 0;

  arma::uvec medoid(k); medoid.zeros();

  bool index_temp_med = true;
  bool value_temp_med = true;

  double ind = 0;
  double ind_max = 0;

  int nan_check = 1;

  while(prova && (iter < maxit))
  {

    iter++;
    U_old = U;
    medoid.zeros();

    for(int c = 0; c < k; c++)
    {
      min_med_old = min_med_const;

      for(int i=0; i<n; i++)
      {
        nan_check = arma::is_finite(U.row(i));
        if(nan_check == 0){
          break;
        }
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

      H.row(c) = data.row(medoid(c));

    }

    D = euclidean_distance_medoid(data, H, n, k);

    U = memb_degree_medoid(D,medoid, m, n, k, p);

    prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

  }

  if(nan_check == 0)
  {
    ind = arma::datum::nan;
    warning("Some membership degrees are NaN (Suggestion: Increase the number of starting points RS)");

  }else{

    value = accu(pow(U,m)%D);
    ind = indices(index, data, U, H, m,  n, k, p,  exp(1.0), alpha);

    // if(index == "PE" || index == "PC" || index == "XB")
    //   {ind_max = -ind;
    //   }else{
    //     ind_max = ind;
    //   }

    if(index == "PE" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }
  }

  return List::create(Rcpp::Named("H") = H,
                      Rcpp::Named("U") = U,
                      Rcpp::Named("D") = D,
                      Rcpp::Named("medoid") = medoid + 1,
                      Rcpp::Named("iter") = iter,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k);

}




/////////////////////////////////////////////////////////


///////////////////////////////////
//        NOISE MEDOID           //
///////////////////////////////////

// [[Rcpp::export]]

List mainFKM_med_noise(arma::mat data,
                       double m,
                       double delta,
                       unsigned int n,
                       unsigned int p,
                       unsigned int k,
                       unsigned rs,
                       double conv,
                       unsigned int maxit,
                       std::string index,
                       double alpha) {

  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U(n,k); U.zeros();

  double func = 0;
  double func_opt = 0;

  arma::mat H_opt(k, p); H.zeros();
  arma::mat U_opt(n,k); U.zeros();
  arma::mat U_old(n,k); U_old.zeros();
  arma::uvec medoid_opt(k); medoid_opt.zeros();
  arma::uword mm = 0;

  double min_med_const = pow(10.0,5.0) * accu(pow(data,2.0));
  double min_med_old = 0;
  double min_med = 0;

  arma::uvec medoid(k); medoid.zeros();

  bool index_temp_med = true;
  bool value_temp_med = true;

  arma::rowvec d(n); d.zeros();

  bool index_temp_med_U = true;
  arma::uvec find_temp(1); find_temp.zeros();

  arma::vec Uout(n); Uout.zeros();
  double q1 = 0;
  double q2 = 0;

  double deltasq = pow(delta,2.0);

  double ind = 0;
  double ind_max = 0;

  int nan_check = 1;

  for(int r=0; r<rs; r++)
  {

    int iter = 0;
    U = unifInit(n,k);
    U_old = U;

    bool prova = true;

    while(prova && (iter < maxit))
    {

      iter++;
      U_old = U;

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

        H.row(c) = data.row(medoid(c));

      }

      D = euclidean_distance_medoid(data, H, n, k);

      for(int i=0; i<n;i++)
      {
        nan_check = arma::is_finite(U.row(i));
        if(nan_check == 0){
          break;
        }

        d = D.row(i);
        index_temp_med_U = Match(i,medoid);

        if(index_temp_med_U)
        {
          if(d.min() == 0)
          {

            mm = d.index_min();
            U(i,mm) = 1;

          }else{

            for(int j=0; j<k;j++)
            {

              q1 = pow(1/arma::as_scalar(D(i,j)),1/(m-1.0)) / sum(pow(1/D.row(i),1/(m-1.0)));
              q2 = pow(arma::as_scalar(D(i,j)),1/(m-1.0)) / pow(delta,2/(m-1.0));
              U(i,j) = 1/((1/q1)+q2);
            }
          }
        }else{


          find_temp = arma::find(medoid == i);
          U.row(i).zeros();
          U(i,find_temp(0)) = 1;
        }

        Uout(i) = 1 - sum(U.row(i));
        Uout(i) = pow(Uout(i),m) * deltasq;
      }

      prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

    }

    if(nan_check == 0)
    {
      func = arma::datum::nan;
    }else{

      func = accu(pow(U,m)%D) + sum(Uout);
    }

    it(r) = iter;
    value(r) = func;

    if(arma::is_finite(func) == true)
    {
      if ((r == 0) | (func < func_opt))
      {
        U_opt=U;
        H_opt=H;
        medoid_opt = medoid;
        func_opt=func;
      }
    }
  }

  if(arma::is_finite(func_opt)){

    ind = indices(index, data, U_opt, H_opt, m,  n, k, p,  exp(1.0), alpha);

    // if(index == "PE" || index == "PC" || index == "XB")
    //   {ind_max = -ind;
    //   }else{
    //     ind_max = ind;
    //   }

    if(index == "PE" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }
  }else{

    ind = arma::datum::nan;
    warning("Some membership degrees are NaN (Suggestion: Increase the number of starting points RS)");

  }


  return List::create(Rcpp::Named("H") = H_opt,
                      Rcpp::Named("U") = U_opt,
                      Rcpp::Named("medoid") = medoid_opt + 1,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k);

}

// [[Rcpp::export]]

List mainFKM_med_noise_U(arma::mat data,
                         double m,
                         double delta,
                         unsigned int n,
                         unsigned int p,
                         unsigned int k,
                         arma::mat U,
                         double conv,
                         unsigned int maxit,
                         std::string index,
                         double alpha) {

  int iter = 0;
  arma::mat H(k, p); H.zeros();
  arma::mat D(n, k); D.zeros();
  arma::mat U_old = U;

  bool prova = true;

  double value = 0;

  double min_med_const = pow(10.0,5.0) * accu(pow(data,2.0));
  double min_med_old = 0;
  double min_med = 0;

  arma::uvec medoid(k); medoid.zeros();

  bool index_temp_med = true;
  bool value_temp_med = true;

  arma::rowvec d(n); d.zeros();

  bool index_temp_med_U = true;
  arma::uvec find_temp(1); find_temp.zeros();
  arma::uword mm = 0;

  arma::vec Uout(n); Uout.zeros();
  double q1 = 0;
  double q2 = 0;

  double deltasq = pow(delta,2.0);

  double ind = 0;
  double ind_max = 0;

  int nan_check = 1;

  while(prova && (iter < maxit))
  {

    iter++;
    U_old = U;

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

      H.row(c) = data.row(medoid(c));

    }

    D = euclidean_distance_medoid(data, H, n, k);

    for(int i=0; i<n;i++)
    {

      nan_check = arma::is_finite(U.row(i));
      if(nan_check == 0){
        break;
      }
      d = D.row(i);
      index_temp_med_U = Match(i,medoid);

      if(index_temp_med_U)
      {
        if(d.min() == 0)
        {

          mm = d.index_min();
          U(i,mm) = 1;

        }else{

          for(int j=0; j<k;j++)
          {

            q1 = pow(1/arma::as_scalar(D(i,j)),1/(m-1.0)) / sum(pow(1/D.row(i),1/(m-1.0)));
            q2 = pow(arma::as_scalar(D(i,j)),1/(m-1.0)) / pow(delta,2/(m-1.0));
            U(i,j) = 1/((1/q1)+q2);
          }
        }
      }else{


        find_temp = arma::find(medoid == i);
        U.row(i).zeros();
        U(i,find_temp(0)) = 1.0;

      }

      Uout(i) = 1.0 - sum(U.row(i));
      Uout(i) = pow(Uout(i),m) * deltasq;
    }

    prova = arma::as_scalar(accu(abs(U_old - U)))>conv;

  }

  if(nan_check == 0)
  {
    ind = arma::datum::nan;
    warning("Some membership degrees are NaN (Suggestion: Increase the number of starting points RS)");


  }else{

    value = accu(pow(U,m)%D) + sum(Uout);
    ind = indices(index, data, U, H, m,  n, k, p,  exp(1.0), alpha);

    // if(index == "PE" || index == "PC" || index == "XB")
    //   {ind_max = -ind;
    //   }else{
    //     ind_max = ind;
    //   }

    if(index == "PE" || index == "XB")
    {ind_max = -ind;
    }else{
      ind_max = ind;
    }
  }

  return List::create(Rcpp::Named("H") = H,
                      Rcpp::Named("U") = U,
                      Rcpp::Named("medoid") = medoid + 1,
                      Rcpp::Named("iter") = iter,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k);

}



//////////////////////////

///////////////////////////////////
//        NEFRC                   //
///////////////////////////////////


// [[Rcpp::export]]

List mainnefrc(arma::mat D,
               double m,
               unsigned int n,
               unsigned int k,
               unsigned int rs,
               double conv,
               unsigned int maxit,
                 std::string index,
                 double alpha) {

  int iter = 0;
  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat U(n,k); U.zeros();
  arma::mat U_old = U;
  arma::mat A(n,k); A.zeros();
  arma::mat b(n,k); b.zeros();
  arma::rowvec B1(k); B1.zeros();
  arma::rowvec Unc(k); Unc.zeros();
  arma::vec Unc1(k); Unc1.zeros();
  double mv = 0;
  arma::uvec mp; mp.zeros();
  double dens = 0;
  arma::rowvec u(k); u.zeros();

  double ind = 0;
  double ind_max = 0;

  double func = 0;
  double func_opt = 0;
  double func_old = 0;

  double first_temp = 0;
  double second_temp = 0;

  arma::mat U_opt(n,k); U_opt.zeros();

  bool prova = true;


  for(int r=0; r<rs; r++)
  {
    iter = 0;
    U = unifInit(n,k);

    prova = true;

    for(int i=0; i<n;i++){
      for(int  j=0; j<k;j++)
      {

        first_temp  = m * arma::as_scalar((pow(U.col(j),m).t() * D.col(i))) / arma::as_scalar(sum(pow(U.col(j),m)));
        second_temp = m * accu((pow(U.col(j),m) * pow(U.col(j),m).t()) % D) / (2 * pow(arma::as_scalar(sum(pow(U.col(j),m))),2.0));
        A(i,j) = first_temp - second_temp;
        b(i,j) = A(i,j) * pow(U(i,j),(m-2.0));
      }
    }

    while(prova && (iter < maxit))
    {
      iter++;
      U_old = U;
      func_old = func;


      for(int i=0; i<n;i++){
        for(int  j=0; j<k;j++)
        {

          first_temp  = m * arma::as_scalar((pow(U.col(j),m).t() * D.col(i))) / arma::as_scalar(sum(pow(U.col(j),m)));
          second_temp = m * accu((pow(U.col(j),m) * pow(U.col(j),m).t()) % D) / (2 * pow(arma::as_scalar(sum(pow(U.col(j),m))),2.0));
          A(i,j) = first_temp - second_temp;
          b(i,j) = A(i,j) * pow(U(i,j),(m-2.0));
        }

        if((double)std::abs(arma::min(b.row(i))) > pow(10.0,-9.0))
        {
          B1 = 1/b.row(i);
          dens = sum(B1.elem(arma::find(B1 > 0)));
          for(int j=0; j < k; j++)
          {
            U(i,j) = (1/arma::as_scalar(b(i,j))) / dens;
            if(U(i,j) < 0)
            {
              U(i,j) = 0;
            }

          }
        }else{

          mv = arma::min(b.row(i));
          mp = arma::find(b.row(i) == mv);
          arma::rowvec uu(k); uu.zeros();
          uu.elem(mp).ones();
          U.row(i) = uu;

        }
        // if(!(is_finite(imag(arma::conv_to<arma::cx_vec>::from(U(i,j))))) || !is_finite(arma::imag(arma::conv_to<arma::cx_vec>::from(U(i,j)))))
        // {
        //   U(i,j) = 0;
        // }
      }
      func = 0;

      for(int j = 0; j<k;j++)
      {
        func += accu((pow(U.col(j),m) * pow(U.col(j),m).t())%D)/ (2*arma::as_scalar(sum(pow(U.col(j),m))));;
      }

      if(!arma::is_finite(func))
      {
        U = U_old;
        func = func_old;
      }

      prova = arma::as_scalar(accu(pow(U_old - U,2.0)))>conv;
    }
    it(r) = iter;
    value(r) = func;

    if(arma::is_finite(func) == true)
    {
      if ((r == 0) | (func < func_opt))
      {
        U_opt=U;
        func_opt=func;
      }
    }
  }

  arma::mat H(1,1); H.ones();
  ind = indices(index, D, U_opt, H, m,  n, k, n,  exp(1.0), alpha, true);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }

  return List::create(Rcpp::Named("U") = U_opt,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max,
                      Rcpp::Named("k") = k);

}


// [[Rcpp::export]]

List mainnefrc_U(arma::mat D,
               arma::mat U,
               double m,
               unsigned int n,
               unsigned int k,
               double conv,
               unsigned int maxit,
               std::string index,
               double alpha) {

  int iter = 0;
  arma::vec value; value.zeros();
  arma::vec it; it.zeros();
  arma::mat U_old = U;
  arma::mat A(n,k); A.zeros();
  arma::mat b(n,k); b.zeros();
  //arma::vec bnc(n); bnc.zeros();
  arma::rowvec B1(k); B1.zeros();
  // arma::rowvec Unc(k); Unc.zeros();
  // arma::vec Unc1(k); Unc1.zeros();
  double mv = 0;
  arma::uvec mp; mp.zeros();
  double dens = 0;
  arma::rowvec u(k); u.zeros();

  double ind = 0;
  double ind_max = 0;

  double func = 0;
  double func_old = 0;

  double first_temp = 0;
  double second_temp = 0;


  bool prova = true;


    while(prova && (iter < maxit))
    {
      iter++;
      U_old = U;
      func_old = func;


      for(int i=0; i<n;i++){
        for(int  j=0; j<k;j++)
        {

          first_temp  = m * arma::as_scalar((pow(U.col(j),m).t() * D.col(i))) / arma::as_scalar(sum(pow(U.col(j),m)));
          second_temp = m * accu((pow(U.col(j),m) * pow(U.col(j),m).t()) % D) / (2 * pow(arma::as_scalar(sum(pow(U.col(j),m))),2.0));
          A(i,j) = first_temp - second_temp;
          b(i,j) = A(i,j) * pow(U(i,j),(m-2));
        }

        if((double)std::abs(arma::min(b.row(i))) > pow(10.0,-9.0))
        {
          B1 = 1/b.row(i);
          dens = sum(B1.elem(arma::find(B1 > 0)));
          for(int j=0; j < k; j++)
          {
            U(i,j) = (1/arma::as_scalar(b(i,j))) / dens;
            if(U(i,j) < 0)
            {
              U(i,j) = 0;
            }

          }
        }else{

          mv = arma::min(b.row(i));
          mp = arma::find(b.row(i) == mv);
          arma::rowvec uu(k); uu.zeros();
          uu.elem(mp).ones();
          U.row(i) = uu;

        }

        // if(!(is_finite(imag(arma::conv_to<arma::cx_vec>::from(U(i,j))))) || !is_finite(arma::imag(arma::conv_to<arma::cx_vec>::from(U(i,j)))))
        // {
        //   U(i,j) = 0;
        // }
      }
      func = 0;

      for(int j = 0; j<k;j++)
      {
        func += accu((pow(U.col(j),m) * pow(U.col(j),m).t())%D)/ (2*arma::as_scalar(sum(pow(U.col(j),m))));;
      }

      if(!arma::is_finite(func))
      {
        U = U_old;
        func = func_old;
      }

      prova = arma::as_scalar(accu(pow(U_old - U,2.0)))>conv;
    }


    it = iter;
    value = func;


  arma::mat H(1,1); H.ones();
  ind = indices(index, D, U, H, m,  n, k, n,  exp(1.0), alpha, true);
  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }

  return List::create(Rcpp::Named("U") = U,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max);

}

//////////////////////////

///////////////////////////////////
//        RNEFRC                 //
///////////////////////////////////

// [[Rcpp::export]]

List mainrnefrc(arma::mat D,
                double m,
                double delta,
                unsigned int n,
                unsigned int k,
                unsigned int rs,
                double conv,
                unsigned int maxit,
                  std::string index,
                  double alpha) {

  int iter = 0;
  arma::vec value(rs); value.zeros();
  arma::vec it(rs); it.zeros();
  arma::mat U(n,k); U.zeros();
  arma::mat U_old = U;
  arma::mat A(n,k); A.zeros();
  arma::mat b(n,k); b.zeros();
  arma::rowvec B1(k); B1.zeros();
  arma::vec Unc(k); Unc.zeros();
  double B2 = 0;
  double mv = 0;
  arma::uvec mp; mp.zeros();
  double dens = 0;
  arma::rowvec u(k); u.zeros();
  arma::vec bnc(n); bnc.zeros();


  double ind = 0;
  double ind_max = 0;

  double func = 0;
  double func_opt = 0;
  double func_old = 0;

  double first_temp = 0;
  double second_temp = 0;

  arma::mat U_opt(n,k); U_opt.zeros();

  bool prova = true;


  for(int r=0; r<rs; r++)
  {
    iter = 0;
    U = unifInit(n,k);

    prova = true;

    while(prova && (iter < maxit))
    {
      iter++;
      U_old = U;
      func_old = func;


      for(int i=0; i<n;i++){
        for(int  j=0; j<k;j++)
        {

          first_temp  = m * arma::as_scalar((pow(U.col(j),m).t() * D.col(i))) / arma::as_scalar(sum(pow(U.col(j),m)));
          second_temp = m * accu((pow(U.col(j),m) * pow(U.col(j),m).t()) % D) / (2 * pow(arma::as_scalar(sum(pow(U.col(j),m))),2.0));
          A(i,j) = first_temp - second_temp;
          b(i,j) = A(i,j) * pow(U(i,j),(m-2.0));
        }

        bnc(i) = m * delta / 2 * pow(sum(U.row(i)),(m-2.0));


        if((double)std::abs(arma::min(b.row(i))) > pow(10.0,-9.0))
        {
          B1 = 1/b.row(i);
          B2 = 1/bnc(i);

          dens = sum(B1.elem(arma::find(B1 > 0)));
          if(B2 > 0)
          {
            dens  += B2;
          }
          for(int j=0; j < k; j++)
          {
            U(i,j) = (1/arma::as_scalar(b(i,j))) / dens;
            if(U(i,j) < 0)
            {
              U(i,j) = 0;
            }

          }
        }else{

          mv = arma::min(b.row(i));
          mp = arma::find(b.row(i) == mv);
          arma::rowvec uu(k); uu.zeros();
          uu.elem(mp).ones();
          U.row(i) = uu;

        }

        // if(!(is_finite(imag(arma::conv_to<arma::cx_vec>::from(U(i,j))))) || !is_finite(arma::imag(arma::conv_to<arma::cx_vec>::from(U(i,j)))))
        // {
        //   U(i,j) = 0;
        // }
      }

      Unc = 1 - sum(U,1);


      func = 0;

      for(int j = 0; j<k;j++)
      {
        func += accu((pow(U.col(j),m) * pow(U.col(j),m).t())%D)/ (2*arma::as_scalar(sum(pow(U.col(j),m))));;
      }

      func += accu((pow(Unc,m) * pow(Unc,m).t()) * delta) / (2*arma::as_scalar(sum(pow(Unc,m))));
      if(!arma::is_finite(func))
      {
        U = U_old;
        func = func_old;
      }

      prova = arma::as_scalar(accu(pow(U_old - U,2.0)))>conv;
    }


    it(r) = iter;
    value(r) = func;

    if(arma::is_finite(func) == true)
    {
      if ((r == 0) | (func < func_opt))
      {
        U_opt=U;
        func_opt=func;
      }
    }
  }

  arma::mat H(1,1); H.ones();
  ind = indices(index, D, U_opt, H, m,  n, k, n,  exp(1.0), alpha,true);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }

  return List::create(Rcpp::Named("U") = U_opt,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max);

}

// [[Rcpp::export]]

List mainrnefrc_U(arma::mat D,
                  arma::mat U,
                  double m,
                  double delta,
                  unsigned int n,
                  unsigned int k,
                  double conv,
                  unsigned int maxit,
                  std::string index,
                  double alpha) {

  int iter = 0;
  arma::vec value; value.zeros();
  arma::vec it; it.zeros();
  arma::mat U_old = U;
  arma::mat A(n,k); A.zeros();
  arma::mat b(n,k); b.zeros();
  arma::rowvec B1(k); B1.zeros();
  arma::vec Unc(k); Unc.zeros();
  double B2 = 0;
  double mv = 0;
  arma::uvec mp; mp.zeros();
  double dens = 0;
  arma::rowvec u(k); u.zeros();
  arma::vec bnc(n); bnc.zeros();


  double ind = 0;
  double ind_max = 0;

  double func = 0;
  double func_old = 0;

  double first_temp = 0;
  double second_temp = 0;

  // arma::mat U_opt(n,k); U_opt.zeros();

  bool prova = true;


  while(prova && (iter < maxit))
  {
    iter++;
    U_old = U;
    func_old = func;


    for(int i=0; i<n;i++){
      for(int  j=0; j<k;j++)
      {

        first_temp  = m * arma::as_scalar((pow(U.col(j),m).t() * D.col(i))) / arma::as_scalar(sum(pow(U.col(j),m)));
        second_temp = m * accu((pow(U.col(j),m) * pow(U.col(j),m).t()) % D) / (2 * pow(arma::as_scalar(sum(pow(U.col(j),m))),2.0));
        A(i,j) = first_temp - second_temp;
        b(i,j) = A(i,j) * pow(U(i,j),(m-2.0));
      }

      bnc(i) = m * delta / 2 * pow(sum(U.row(i)),(m-2.0));


      if((double)std::abs(arma::min(b.row(i))) > pow(10.0,-9.0))
      {
        B1 = 1/b.row(i);
        B2 = 1/bnc(i);

        dens = sum(B1.elem(arma::find(B1 > 0)));
        if(B2 > 0)
        {
          dens  += B2;
        }
        for(int j=0; j < k; j++)
        {
          U(i,j) = (1/arma::as_scalar(b(i,j))) / dens;
          if(U(i,j) < 0)
          {
            U(i,j) = 0;
          }

        }
      }else{

        mv = arma::min(b.row(i));
        mp = arma::find(b.row(i) == mv);
        arma::rowvec uu(k); uu.zeros();
        uu.elem(mp).ones();
        U.row(i) = uu;

      }

      // if(!(is_finite(imag(arma::conv_to<arma::cx_vec>::from(U(i,j))))) || !is_finite(arma::imag(arma::conv_to<arma::cx_vec>::from(U(i,j)))))
      // {
      //   U(i,j) = 0;
      // }
    }

    Unc = 1 - sum(U,1);


    func = 0;

    for(int j = 0; j<k;j++)
    {
      func += accu((pow(U.col(j),m) * pow(U.col(j),m).t())%D)/ (2*arma::as_scalar(sum(pow(U.col(j),m))));;
    }

    func += accu((pow(Unc,m) * pow(Unc,m).t()) * delta) / (2*arma::as_scalar(sum(pow(Unc,m))));
    if(!arma::is_finite(func))
    {
      U = U_old;
      func = func_old;
    }

    prova = arma::as_scalar(accu(pow(U_old - U,2.0)))>conv;
  }


  it = iter;
  value = func;


  arma::mat H(1,1); H.ones();
  ind = indices(index, D, U, H, m,  n, k, n,  exp(1.0), alpha);

  // if(index == "PE" || index == "PC" || index == "XB")
  //   {ind_max = -ind;
  //   }else{
  //     ind_max = ind;
  //   }

  if(index == "PE" || index == "XB")
  {ind_max = -ind;
  }else{
    ind_max = ind;
  }

  return List::create(Rcpp::Named("U") = U,
                      Rcpp::Named("iter") = it,
                      Rcpp::Named("value") = value,
                      Rcpp::Named("k") = k,
                      Rcpp::Named("index") = ind,
                      Rcpp::Named("index_max") = ind_max);

}
