#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::depends(RcppArmadillo)]]

//'Fast C++ computation of the final posterior probabilities in the E-M Winkler's method
//'
//'@param agreemat binary sparse matrix of dimensions \code{N x K} containing the agreement rows for each pair of potential matches.
//'@param m vector of length \code{K} containing the agreement weights.
//'@param u vector of length \code{K} containing the disagreement weights.
//'@param nA integer indicating the number of observations to be matched.
//'@param nB integer indicating the number of observations to be matched with.
//'
//'@export
// [[Rcpp::export]]
arma::mat matchingScore_C(arma::mat agreemat, arma::vec m, arma::vec u, int nA, int nB){
  
  int K = agreemat.n_cols;  
  int N = agreemat.n_rows;
  
  //vec temp_res = vec(N);
  mat res(nA, nB);
  
  //umat agreemat_neg = 1 - agreemat;
  
  for(int j=0; j<nB; j++){
    for(int i=0; i<nA; i++){
      int l = j*nA + i;
      res(i,j) = 0.0;
      for(int k=0; k<K; k++){
        //temp_res(j) = sum(pow(log(m)-log(u), round(agreemat.row(j))) * pow(log(1-m)-log(1-u), 1-round(agreemat.row(j))));
        int expo(agreemat(l,k));
        res(i,j) += pow(log(m(k))-log(u(k)), expo) * pow(log(1-m(k))-log(1-u(k)), 1-expo);
      }
    }
  }
  
  //int j1 = i*nA;
  //int j2 = i*nA + nA-1;
  //res.col(i) = temp_res.subvec(j1,j2);
  
  return(res);
}


//'@rdname matchingScore_C
//'@param mat_A a \code{nB x K} matrix of the observations to be matched.
//'@param mat_B a \code{nA x K} matrix of the database into which a match is looked for.
//'@description \code{matchingScore_C_sparse_big} implements a version using sparse matrices. It has a better 
//'management of memory but is a little bit slower (indicated for big matrices)
//'@export
// [[Rcpp::export]]
arma::mat matchingScore_C_sparse_big(arma::mat mat_A, arma::mat mat_B, arma::vec m, arma::vec u){
  
  int K = mat_A.n_cols;  
  int nA = mat_A.n_rows;
  int nB = mat_B.n_rows;
  int N = nA*nB;
  
  //vec temp_res = vec(N);
  mat res(nA, nB);
  
  //umat agreemat_neg = 1 - agreemat;
  
  for(int j=0; j<nB; j++){
    for(int i=0; i<nA; i++){
      int l = j*nA + i;
      res(i,j) = 0.0;
      for(int k=0; k<K; k++){
        //temp_res(j) = sum(pow(log(m)-log(u), round(agreemat.row(j))) * pow(log(1-m)-log(1-u), 1-round(agreemat.row(j))));
        int agree_neg(as_scalar(mat_A(i,k) - mat_B(j,k)));
        int expo(1 - abs(agree_neg));
        res(i,j) += pow(log(m(k))-log(u(k)), expo) * pow(log(1.0-m(k))-log(1.0-u(k)), 1-expo);
      }
    }
  }
  
  //int j1 = i*nA;
  //int j2 = i*nA + nA-1;
  //res.col(i) = temp_res.subvec(j1,j2);
  
  return(res);
}
