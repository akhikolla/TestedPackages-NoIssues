// [[Rcpp::depends(bigstatsr, rmio, RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <bigstatsr/BMAcc.h>
#include <iostream>
#include <stdio.h>
#include <strings.h>

using namespace arma;
using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

//' @title Kernel parameter estimation
//' @description Kernel parameter estimation by averaging the distances to the closest neighbors.
//' @param X A Filebacked Big Matrix \eqn{n \times N}.
//' @param size Neighborhood size.
//' @param kernel_type Kernel function type. Available types are c("Gaussian", "Laplacian").
//' @return The estimated kernel parameter.
//' @export
// [[Rcpp::export]]
double gamma_estimation(Environment X, int size, const char* kernel_type){
  XPtr<FBM_RW> xp_X = X["address_rw"];
  BMAcc_RW<double> mat_X(xp_X);
  size_t N = mat_X.nrow();
  size_t n = mat_X.ncol();
  double rowNeighborhood = 0;
  NumericVector Neighbors(n);
  NumericVector out(size);
  for(size_t i=0; i< N; i++){
    for(size_t j=0; j< n; j++){
      Neighbors[j] = mat_X(i,j);
    }
    std::sort(Neighbors.begin(), Neighbors.end());
    std::copy(Neighbors.begin(), Neighbors.begin() + size, out.begin());
    
  //  if(std::strcmp(kernel_type, "Gaussian")==0){
      rowNeighborhood += sum(out*out);  
    // }else{
    //    rowNeighborhood += sum(out);
    // }
    // 
  }
  if(std::strcmp(kernel_type, "Gaussian")==0){
    return (N*size)/(2.0*rowNeighborhood);
  }else{
    return (std::sqrt(N*size))/(std::sqrt(rowNeighborhood));
  }
}

struct CumsumParallel : public Worker {
  BMAcc_RW<double> M_input;
  BMAcc_RW<double> M_output;
  size_t n;
  CumsumParallel(BMAcc_RW<double> M_input, BMAcc_RW<double> M_output)
    : M_input(M_input),M_output(M_output) {
    n = M_input.nrow();
    }
  void operator()(size_t begin, size_t end) {
    for(size_t j = begin; j< end; j++) {
      double sum = 0;
      for (size_t i = 0; i < n; i++) {
        sum += M_input(i,j);
        M_output(i,j) = sum;
      }
    }
  }
};
//' @title Cumulative sum computation
//' @description Parallel implementation of the cumulative sum of the matrix columns.
//' @param X A Filebacked Big Matrix n x N.
//' @param A_cumsum A Filebacked Big Matrix n x N, where cummulative sums are stored.
//' @return NULL
//' @export
// [[Rcpp::export]]
void cumsum_parallel(Environment X, Environment A_cumsum){
  XPtr<FBM_RW> xp_M_in= X["address_rw"];
  BMAcc_RW<double> M_input(xp_M_in);

  XPtr<FBM_RW> xp_M_out= A_cumsum["address_rw"];
  BMAcc_RW<double> M_output(xp_M_out);

  CumsumParallel cumsum_rcpp(M_input, M_output);
  parallelFor(0, M_input.ncol(), cumsum_rcpp);
}

struct Wasserstein_Parallel : public Worker {
  BMAcc_RW<double> M_input;
  BMAcc_RW<double> M_output;
  IntegerVector& set_c;
  size_t c, n;
  Wasserstein_Parallel(BMAcc_RW<double> M_input, BMAcc_RW<double> M_output, IntegerVector& set_c)
    : M_input(M_input),M_output(M_output), set_c(set_c) {
    c = set_c.size();
    n = M_input.nrow();
  }
  void operator()(size_t begin, size_t end) {
    for(size_t i = begin; i< end; i++) {

      for(size_t j = 0; j < c; j++){

        double dist = 0;
        size_t row_current = set_c[j] -1;
        double last_current = M_input(n-1, row_current);
        double last =  M_input( n-1, i);

        for( size_t k =0; k < n; k++){
          dist += std::abs(M_input(k, row_current)/last_current - M_input(k, i)/last);
        }
        M_output(i,j) = dist;
      }
    }
  }
};
//' @title Wasserstein-1 distance 
//' @param X A Filebacked Big Matrix n x N.
//' @param C A Filebacked Big Matrix N x l, which stores the Wasserstein distances. 
//' @param set_c Column index vector. The data vector indices for which the Wasserstein distances are computed. 
//' @return NULL
//' @details The Wasserstein-1 distances are computed between the data vectors from \code{set_c} and all columns of X. 
//' @export
// [[Rcpp::export]]
void W1_parallel(Environment X, Environment C, IntegerVector set_c){
  XPtr<FBM_RW> xp_M_in= X["address_rw"];
  BMAcc_RW<double> M_input(xp_M_in);
  
  XPtr<FBM_RW> xp_M_out= C["address_rw"];
  BMAcc_RW<double> M_output(xp_M_out);

  Wasserstein_Parallel W1_distance(M_input, M_output, set_c);
  parallelFor(0, M_input.ncol(), W1_distance);
}

struct Euclidean_Parallel : public Worker {
  BMAcc_RW<double> M_input;
  BMAcc_RW<double> M_output;
  IntegerVector& set_c;
  size_t c, n;
  Euclidean_Parallel(BMAcc_RW<double> M_input, BMAcc_RW<double> M_output, IntegerVector& set_c)
    : M_input(M_input),M_output(M_output), set_c(set_c) {
    c = set_c.size();
    n = M_input.nrow();
  }
  void operator()(size_t begin, size_t end) {
    for(size_t i = begin; i< end; i++) {
      
      for(size_t j = 0; j < c; j++){
        
        double dist = 0;
        size_t row_current = set_c[j] -1;
       // double last_current = M_input(n-1, row_current);
        //double last =  M_input( n-1, i);
        
        for( size_t k =0; k < n; k++){
          dist += std::pow(M_input(k, row_current) - M_input(k, i), 2.);
        }
        M_output(i,j) = std::sqrt(dist);
      }
    }
  }
};
//' @title Euclidean distance 
//' @param X A Filebacked Big Matrix \eqn{n \times N}.
//' @param C A Filebacked Big Matrix \eqn{N \times l}, which stores the Euclidean distances. 
//' @param set_c Column index vector. The data vector indices for which the Euclidean distances are computed. 
//' @return NULL
//' @details The Euclidean distances are computed between the data vectors from \code{set_c} and all columns of X. 
//' @export
// [[Rcpp::export]]
void E_parallel(Environment X, Environment C, IntegerVector set_c){
  XPtr<FBM_RW> xp_M_in= X["address_rw"];
  BMAcc_RW<double> M_input(xp_M_in);
  
  XPtr<FBM_RW> xp_M_out= C["address_rw"];
  BMAcc_RW<double> M_output(xp_M_out);
  
  Euclidean_Parallel E_distance(M_input, M_output, set_c);
  parallelFor(0, M_input.ncol(), E_distance);
}
