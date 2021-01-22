#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


//To remember PKG_CXXFLAGS = -DARMA_64BIT_WORD

//[[Rcpp::depends(RcppArmadillo)]]



//'Fast C++ implementation of agreement vector for the element-wise comparison of 2 matrices
//'
//'@param mat_A a \code{nB x K} matrix of the observations to be matched. Must be integers.
//'
//'@param mat_B a \code{nA x K} matrix of the database into which a match is looked for. Must be integers.
//'
//'@aliases agree_C_sparse agree_C
//'
//'@examples
//'mat1 <- matrix(round(rnorm(n=1000, sd=1.2)), ncol=10, nrow=100)
//'mat2 <- rbind(mat1[1:10, ],
//'              matrix(round(rnorm(n=900, sd=1.2)), ncol=10, nrow=90)
//'              )
//'rownames(mat1) <- paste0("A", 1:nrow(mat1))
//'rownames(mat1) <- paste0("B", 1:nrow(mat1))
//'mat1 <- 1*(mat1>1)
//'mat2 <- 1*(mat2>1)
//'
//'@export
// [[Rcpp::export]]
arma::mat agree_C(arma::mat mat_A, 
                  arma::mat mat_B){
  
  int K = mat_A.n_cols;  
  int nA = mat_A.n_rows;
  int nB = mat_B.n_rows;
  int N = nA*nB;
  
  mat res(N,K);
  
  for(int j=0; j<nB; j++){
    int j1 = j*nA;
    int j2 = j*nA + nA-1;
    res.rows(j1,j2) = 1.0 - abs(mat_A.each_row() - mat_B.row(j));
  }
  
  return(res);
}


//'@rdname agree_C
//'
//'@description \code{agree_C_sparse} uses sparse matrices.
//'
//'@export
// [[Rcpp::export]]
arma::sp_mat agree_C_sparse(arma::mat mat_A, 
                            arma::mat mat_B){
  
  int K = mat_A.n_cols;  
  int nA = mat_A.n_rows;
  int nB = mat_B.n_rows;
  int N = nA*nB;
  
  sp_mat res(N,K);
  
  for(int j=0; j<nB; j++){
    int j1 = j*nA;
    int j2 = j*nA + nA-1;
    res.rows(j1,j2) = 1.0 - abs(mat_A.each_row() - mat_B.row(j));
  }
  
  return(res);
}
