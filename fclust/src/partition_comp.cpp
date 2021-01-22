#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


List partition_comp(arma::mat HardClust, arma::mat Fuzzy, std::string t_norm)
{


  arma::mat R = HardClust.t();
  arma::mat Q = Fuzzy.t();

  unsigned int N = R.n_cols;
  unsigned int k = R.n_rows;
  unsigned int v = Q.n_rows;

  arma::mat V(N,N); V.zeros();
  arma::mat X(N,N); X.zeros();
  arma::mat Y(N,N); Y.zeros();
  arma::mat Z(N,N); Z.zeros();

  bool tnorm = true;

  double a = 0;
  double b = 0;
  double c = 0;
  double d = 0;

  double Rand_F = 0;
  double ARand_F = 0;
  double Jaccard_F = 0;

  double m=0;

  if(t_norm == "product") {tnorm = false;}

  for(int j2 = 1; j2<N; j2++)
  {
    for(int j1 = 0; j1<(j2); j1++)
    {

      for(int i = 0; i<k; i++)
      {

        if(tnorm == true){
          m = std::min(R(i,j1),R(i,j2));
          V(j1,j2) = std::max(V(j1,j2),m);
        }
        if(tnorm == false){
          V(j1,j2) = std::max(V(j1,j2),R(i,j1)*R(i,j2));
        }
      }

      for(int i1 = 0; i1<k; i1++)
      {
        for(int i2 = 0; i2<k; i2++)
        {
          if(i1 != i2)
          {
            if(tnorm == true){
              m = std::min(R(i1,j1),R(i2,j2));
              X(j1,j2) = std::max(X(j1,j2),m);
            }
            if(tnorm == false){
              X(j1,j2) = std::max(X(j1,j2),R(i1,j1)*R(i2,j2));
            }
          }
        }
      }

      for(int l = 0; l<v; l++)
      {

        if(tnorm == true){
          m = std::min(Q(l,j1),Q(l,j2));
          Y(j1,j2) = std::max(Y(j1,j2),m);
        }
        if(tnorm == false){
          Y(j1,j2) = std::max(Y(j1,j2),Q(l,j1)*Q(l,j2));
        }
      }

      for(int l1 = 0; l1<v; l1++)
      {
        for(int l2 = 0; l2<v; l2++)
        {
          if(l1 !=l2)
          {
            if(tnorm == true){
              m = std::min(Q(l1,j1),Q(l2,j2));
              Z(j1,j2) = std::max(Z(j1,j2),m);
            }
            if(tnorm == false){
              Z(j1,j2) = std::max(Z(j1,j2),Q(l1,j1)*Q(l2,j2));
            }
          }
        }
      }

      if(tnorm == true)
      {

        a = a + std::min(V(j1,j2),Y(j1,j2));
        b = b + std::min(V(j1,j2),Z(j1,j2));
        c = c + std::min(X(j1,j2),Y(j1,j2));
        d = d + std::min(X(j1,j2),Z(j1,j2));
      }
      if(tnorm == false)
      {

        a = a + (V(j1,j2)*Y(j1,j2));
        b = b + (V(j1,j2)*Z(j1,j2));
        c = c + (X(j1,j2)*Y(j1,j2));
        d = d + (X(j1,j2)*Z(j1,j2));
      }

    }


  }


  Rand_F = (a+d)/(a+b+c+d);
  ARand_F = (2*(a*d-b*c))/((pow(b,2)) + (pow(c,2)) + (2*a*d) + ((a+d) * (c+b)));
  Jaccard_F = a / (a+b+c);

  return List::create(Rcpp::Named("Rand.F") = Rand_F,
                      Rcpp::Named("adjRand.F") = ARand_F,
                      Rcpp::Named("Jaccard.F") = Jaccard_F);

}
