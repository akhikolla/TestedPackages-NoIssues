#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "distance.h"

using namespace Rcpp;



//////////////////////////////////////////////
///   USEFUL FUNCTIONS                       //
//                                          //
//                                          //
//////////////////////////////////////////////


///////////////////////////////////
// REPLACING ELEMENTS IN A VECTOR//
///////////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::rowvec replace(arma::rowvec x, double val, double change) {

  arma::uvec ids = find(x < val);

  x.elem(ids).fill(change);

  return x;
}

////////////////////////////////////////////////////////////
// CHECK IF THE INVERSE CAN BE COMOUTED                   //
// RETURNS A BOLEAN false IF NOT AND THE INVERSE IS true //
//////////////////////////////////////////////////////////

// // [[Rcpp::export]]
//
// arma::mat InvCheck(arma::mat A)
// {
//
//   arma::mat inverse = A;
//
//   try
//   {
//
//     inverse = arma::inv_sympd(A);
//   }
//   catch(...)
//   {
//
//     //A = false;
//   }
//
//   if(is_finite(A) == FALSE)
//   {
//     inverse.reset();
//   }
//   return inverse;
// }

// [[Rcpp::export]]

arma::mat InvCheck(arma::mat A)
{

  arma::mat inverse = A;

  double eps = std::numeric_limits<double>::epsilon();
  double rcondition = 0;

  rcondition = arma::rcond(A);

  if(rcondition > eps)
  {
    inverse = arma::inv(A);
  }else{
    inverse.reset();
  }
  return inverse;
}

////////////////////////////////////////////////////////////
// FINDS ELEMTENS IN A VECTOR, RETURNING A BOOL VALUE    //
//  true if match and false otherwise                   //
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
bool Match(int i, arma::uvec B)
{

  arma::uvec A = arma::find(B == i);

  bool out = A.is_empty();

  return out;
}


////////////////////////////////////////////////////////////
// IF D IS A DISTANCE MATRIX                             //
//////////////////////////////////////////////////////////

// [[Rcpp::export]]

void distCheck(NumericMatrix D, unsigned int n, unsigned int p)
{

  bool diag = true;
  bool sym = true;
  bool pos = true;
  bool finite = true;

  if(n != p)
  {
    stop("D must be a square matrix");
  }

  for(int i=0; i<n;i++)
  {
    for(int j=0; j<i+1; j++)
    {

      finite = (arma::is_finite(D(i,j)) & arma::is_finite(D(j,i)));
      if(finite == false)
      {
        stop("The data set X must not contain NA values");
      }

      sym = (D(i,j) == D(j,i));
      if(sym == false)
      {
        stop("D must be symmetric");
      }

      pos = (D(i,j) >= 0);
      if(pos == false)
      {
        stop("Each value of D must be positive");
      }
    }

    diag = (D(i,i) == 0);
    if(diag == false)
    {
      stop("The diagonal of D must be 0");
    }
  }
}


/////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////
///   CLUSTER VALIDITY INDICES             //
//                                          //
//                                          //
//////////////////////////////////////////////

////////////////////////////////////////////////////////////
// silhouette INDEX FOR  INTERNAL USE                   //
//////////////////////////////////////////////////////////

// [[Rcpp::export]]

double silhouette_internal(arma::mat X, arma::mat U, unsigned int p, unsigned int k, unsigned int n,bool distance = false)
{

  arma::vec memb(n); memb.zeros();
  arma::ivec count(k); count.zeros();

  arma::mat D(n,n); D.zeros();

  arma::vec a(n); a.zeros();
  arma::vec b = a;
  arma::vec sil_obj = a;
  arma::mat B(n,k); B.zeros();

  double sil = 0;

  for(int i = 0; i<n;i++)
  {

    memb(i) = U.row(i).index_max();
    count(memb(i)) += 1;

  }

  if(distance == false)
  {
    for(int i=0;i<(n-1);i++)
    {
      for(int i2=(i+1); i2<n;i2++)
      {
        D(i,i2) = sum(pow(X.row(i)-X.row(i2),2.0));
        D(i2,i) = D(i,i2);

      }
    }
  }else{
    D = X;
  }

  for(int i=0;i<n;i++)
  {
    for(int j=0;j<k;j++)
    {
      for(int i2=0;i2<n;i2++)
      {
        if(memb(i2)==j){B(i,j) += D(i,i2);}
      }
    }
  }
  for(int i = 0; i<n;i++)
  {
    for(int j = 0; j<k;j++)
    {
      if(memb(i) == j)
      {
        if (count(j)!=1)
        {
          B(i,j) = B(i,j)/(count(j)-1.0);
          a(i) = B(i,j);
          B(i,j)=max(B.row(i)) + 1;
        }
      }else{
        B(i,j) = B(i,j)/count(j);
      }
    }

    if(count(memb(i)) != 1)
    {
      b(i)=arma::min(B.row(i));
      sil_obj(i)=(b(i)-a(i))/(std::max(a(i),b(i)));
    }

  }


  sil = arma::mean(sil_obj);
  return sil;
}

////////////////////////////////////////////////////////////
// silhouette INDEX                                      //
//////////////////////////////////////////////////////////

// [[Rcpp::export]]

List silhouette(arma::mat X, arma::mat U, unsigned int p, unsigned int k, unsigned int n, bool distance = false)
{
  arma::uvec memb(n); memb.zeros();
  arma::ivec count(k); count.zeros();

  arma::mat D(n,n); D.zeros();

  arma::vec a(n); a.zeros();
  arma::vec b = a;
  arma::vec sil_obj = a;
  arma::mat B(n,k); B.zeros();

  double sil = 0;

  for(int i = 0; i<n;i++)
  {

    memb(i) = U.row(i).index_max();
    count(memb(i)) += 1;

  }

  if(distance == false)
  {
    for(int i=0;i<(n-1);i++)
    {
      for(int i2=(i+1); i2<n;i2++)
      {
        D(i,i2) = sum(pow(X.row(i)-X.row(i2),2.0));
        D(i2,i) = D(i,i2);

      }
    }

  }else{

    D = X;

  }

  for(int i=0;i<n;i++)
  {
    for(int j=0;j<k;j++)
    {
      for(int i2=0;i2<n;i2++)
      {
        if(memb(i2)==j){B(i,j) += D(i,i2);}
      }
    }
  }
  for(int i = 0; i<n;i++)
  {
    for(int j = 0; j<k;j++)
    {
      if(memb(i) == j)
      {
        if (count(j)!=1)
        {
          B(i,j) = B(i,j)/(count(j)-1);
          a(i) = B(i,j);
          B(i,j)=max(B.row(i)) + 1;
        }
      }else{
        B(i,j) = B(i,j)/count(j);
      }
    }

    if(count(memb(i)) != 1)
    {
      b(i)=arma::min(B.row(i));
      sil_obj(i)=(b(i)-a(i))/(std::max(a(i),b(i)));
    }

  }


  sil = arma::mean(sil_obj);
  return List::create(Rcpp::Named("sil") = sil,
                      Rcpp::Named("sil.obj") = sil_obj);
}


////////////////////////////////////////////////////////////
// FUZZY silhouette INDEX                                //
//////////////////////////////////////////////////////////

// [[Rcpp::export]]

double silhouetteFuzzy(arma::mat X, arma::mat U, double alpha,unsigned int p, unsigned int k, unsigned int n, bool distance = false)
{
  arma::uvec memb(n); memb.zeros();
  arma::ivec count(k); count.zeros();

  arma::mat D(n,n); D.zeros();

  arma::vec a(n); a.zeros();
  arma::vec b = a;
  arma::vec sil_obj = a;
  arma::mat B(n,k); B.zeros();

  double sil = 0;

  arma::vec sil_F(n); sil_F.zeros();
  arma::vec w(n); w.zeros();

  arma::uvec ii; ii.zeros();
  arma::rowvec uu(k-1); uu.zeros();
  arma::rowvec u(k); u.zeros();

  for(int i = 0; i<n;i++)
  {

    memb(i) = U.row(i).index_max();
    count(memb(i)) += 1;

  }

  if(distance == false)
  {
    for(int i=0;i<(n-1);i++)
    {
      for(int i2=(i+1); i2<n;i2++)
      {
        D(i,i2) = sum(pow(X.row(i)-X.row(i2),2.0));
        D(i2,i) = D(i,i2);

      }
    }
  }else{
    D = X;
  }

  for(int i=0;i<n;i++)
  {
    for(int j=0;j<k;j++)
    {
      for(int i2=0;i2<n;i2++)
      {
        if(memb(i2)==j){B(i,j) += D(i,i2);}
      }
    }
  }

  for(int i = 0; i<n;i++)
  {
    for(int j = 0; j<k;j++)
    {
      if(memb(i) == j)
      {
        if (count(j)!=1)
        {
          B(i,j) = B(i,j)/(count(j)-1.0);
          a(i) = B(i,j);
          B(i,j)=max(B.row(i)) + 1;
        }
      }else{
        B(i,j) = B(i,j)/count(j);
      }
    }

    if(count(memb(i)) != 1)
    {
      b(i)=arma::min(B.row(i));
      sil_obj(i)=(b(i)-a(i))/(std::max(a(i),b(i)));
    }

  }

  for(int i = 0; i<n;i++)
  {
    u = arma::sort(U.row(i),"descend");
    uu = u.subvec(1,k-1);
    w(i) = pow(arma::max(u) - arma::max(uu),alpha);

  }

  sil = sum(w.t()*sil_obj)/sum(w);
  return sil;
}

////////////////////////////////////////////////////////////
// PARTITION COEFFICIENT                                   //
//////////////////////////////////////////////////////////
// [[Rcpp::export]]

double partCoef(arma::mat U, unsigned int n)
{

  double out = 0;
  out = accu(pow(U,2.0))/(double)n;
  return out;

}

////////////////////////////////////////////////////////////
// PARTITION ENTROPY INDEX                             //
//////////////////////////////////////////////////////////

// [[Rcpp::export]]

double partEntropy(arma::mat U, double b, unsigned int n)
{

  double out = 0;
  double eps = std::numeric_limits<double>::epsilon();
  arma::mat logBase = U;

  U.elem(find(U < eps)).fill(eps);

  logBase = log(U) / log(b);
  out = -accu(U%logBase)/(double)n;

  return out;

}

////////////////////////////////////////////////////////////
// MODIFIED PARTITION COEFFICIENT                      //
//////////////////////////////////////////////////////////

// [[Rcpp::export]]

double partCoef_mod(arma::mat U, unsigned int n, unsigned int k)
{

  double out = 0;

  out = 1-(double)k/(double)(k-1)*(1-partCoef(U,n));

  return out;

}

////////////////////////////////////////////////////////////
// XIE AND BENI INDEX                                   //
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
double xie_beni(arma::mat X,arma::mat U, arma::mat H, double m, unsigned int n, unsigned int k)
{

  arma::mat D(n,k); D = euclidean_distance(X,H,n,k);
  double distH = pow(10.0,10.0)*accu(pow(H,2.0));

  double out=0;

  for(int i=0; i<(k-1);i++)
  {
    for(int j=(i+1); j<k;j++)
    {
      if (sum(pow(H.row(i)-H.row(j),2.0))<distH)
      {
        distH=sum(pow(H.row(i)-H.row(j),2.0));
      }
    }
  }

  out = accu(pow(U,m)%D)/((double)n*distH);
  return out;

}


////////////////////////////////////////////////////////////
//  INDEX SELECTION                                      //
//////////////////////////////////////////////////////////

// [[Rcpp::export]]
double indices(std::string type, arma::mat X, arma::mat U, arma::mat H, double m, unsigned int n, unsigned int k, unsigned int p, double b, double alpha, bool distance = false)
{

  double value = 0;

  if(type == "PC"){
    value = partCoef(U,n);}
  else if (type == "PE"){
    value = partEntropy(U,b,n);}
  else if (type == "MPC"){
    value = partCoef_mod(U,n,k);}
  else if (type == "SIL"){
    value = silhouette_internal(X,U,p,k,n,distance);}
  else if (type == "SIL.F"){
    value = silhouetteFuzzy(X,U,alpha,p,k,n,distance);}
  else if (type == "XB"){
    value = xie_beni(X,U,H,m,n,k);}
  else {stop("No match names.");}

  return value;
}

//////////////////////////////////////////////
///   RANDOM INITIALIZATION                 //
//GENERATE A INITIAL MEMBERSHIP DEGREE MAT  //
//////////////////////////////////////////////

// // [[Rcpp::export]]
//
// arma::mat unifInit_Rcpp(int n, int d)
// {
//
//   arma::mat unif(n,d); unif.zeros();
//   arma::rowvec unif_temp(d);  unif_temp.zeros();
//   double s = 0;
//
//
//   for(int i=0;i<n;i++)
//   {
//     unif_temp.randu();
//     s = sum(unif_temp);
//     unif.row(i) = unif_temp/s;
//   }
//
//
//   return unif;
// }

// [[Rcpp::export]]

arma::mat unifInit(int n, int d)
{

  arma::mat unif(n,d); unif.zeros();
  arma::colvec unif_temp(n);  unif_temp.zeros();

  for(int i=0;i<d;i++)
  {
    unif_temp.randu();
    unif.col(i) = unif_temp;
  }


  for(int j = 0; j < n; j++)
  {
    unif.row(j) = unif.row(j)/sum(unif.row(j));

  }

  return unif;
}

/////////////////////////////////////////////////////////
