// [[Rcpp::depends(bigstatsr, RcppArmadillo, rmio)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <bigstatsr/BMAcc.h>
#include <iostream>

using namespace arma;
using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

void cumsum_parallel(Environment X, Environment A_cumsum);

struct W1_cl_centr_parallel : public Worker {
  SubBMAcc_RW<double> M_input;
  NumericVector output;
  NumericVector& C_cumsum;
  size_t c, n;
  W1_cl_centr_parallel(SubBMAcc_RW<double> M_input, NumericVector output, NumericVector& C_cumsum)
    : M_input(M_input),output(output), C_cumsum(C_cumsum) {
    c = C_cumsum.size();
    n = M_input.nrow();
  }
  void operator()(size_t begin, size_t end) {
    for(size_t i = begin; i< end; i++) {

        double dist = 0;
        double last_current = M_input(n-1, i);
        double last =  C_cumsum(n-1);

        for( size_t k =0; k < n; k++){
          dist += std::abs(M_input(k, i)/last_current - C_cumsum(k)/last);
        }
        output(i) = dist;
      }
    }
};
//' Parallel computation of the Wasserstein distances between cluster centroid and points in the cluster
//'  
//' @param X_cumsum is a Filebacked data cummulative sum matrix n x N
//' @param rowInd is a vector of the row indeces
//' @param colInd is a vector of the column indeces
//' @param C_cumsum is a vector of the cummulative sum of the cluster centroid
//' @return The vector of the Wasserstein distances
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericVector W1_cl_centr(Environment X_cumsum,const IntegerVector& rowInd, 
                          const IntegerVector& colInd, NumericVector C_cumsum){
  XPtr<FBM_RW> xp_X= X_cumsum["address_rw"];
  SubBMAcc_RW<double> mat_X(xp_X, rowInd -1, colInd-1);
  size_t size_cl = colInd.size();
  NumericVector output(size_cl);

  W1_cl_centr_parallel Dist(mat_X, output, C_cumsum);
  parallelFor(0, size_cl, Dist);
  return output;
}

//' Computation of the Wasserstein distances between cluster centroid and signals in the cluster
//'  
//' @param X_cumsum is a Filebacked data cummulative sum matrix n x N
//' @param rowInd is a vector of the row indeces
//' @param colInd is a vector of the column indeces
//' @param C_cumsum is a vector of the cummulative sum of the cluster centroid
//' @export
//' @return The vector of the Wasserstein distances
//' @keywords internal
// [[Rcpp::export]]
double IntraDist(Environment X_cumsum, const IntegerVector& rowInd, const IntegerVector& colInd, 
                 NumericVector C_cumsum){
  XPtr<FBM_RW> xp_X= X_cumsum["address_rw"];
  SubBMAcc_RW<double> mat_X(xp_X, rowInd -1, colInd-1);
  size_t N = mat_X.ncol();
  size_t n = mat_X.nrow();
  double cl_size = colInd.size();
  double dist=0;
  for(size_t j=0; j< N; j++){
    for(size_t i=0; i< n; i++){
      dist+= abs(mat_X(i,j)/mat_X(n-1,j) - C_cumsum(i)/C_cumsum(n-1)); //pow(mat_X(i,j) - C(i), 2);
    }
  }
  return dist/cl_size;
}
//'@title cumsum_sug
//' @description rcpp cummulative sum of integer vector 
//' @param x is a numeric vector
//' @return A vector of cummulative sum
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericVector cumsum_sug(NumericVector x){
  return cumsum(x);    // compute the result vector and return it
}
//'@title cumsum_Mat
//' @description rcpp cummulative sum of matrix columns 
//' @param C is a numeric matrix
//' @export
//' @return A matrix of cummulative sum
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix cumsum_Mat(NumericMatrix C){
  size_t N = C.ncol();
  size_t n = C.nrow();
  NumericMatrix C_cumsum(n,N);
  for(size_t j=0; j < N; j++){
    C_cumsum(_,j) = cumsum_sug(C(_,j));
  }
  return C_cumsum;
}

//' Computation of the Wasserstein distances between cluster centroids
//'  
//' @param C_cumsum is a Filebacked data cummulative sum matrix of N cluster centroids
//' @export
//' @return N x N matrix of the Wasserstein distances
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix InterDist(NumericMatrix C_cumsum){
  size_t N = C_cumsum.ncol();
  size_t n = C_cumsum.nrow();
  NumericMatrix dist(N,N);
  double sum;

  for(size_t j=0; j< N; j++){
    for(size_t k=0; k<j; k++){
      sum = 0;
      for(size_t i=0; i< n; i++){
        sum += abs(C_cumsum(i,j)/C_cumsum(n-1, j) - C_cumsum(i,k)/C_cumsum(n-1, k)); //
      }
      dist(j,k) = sum;
      dist(k,j) = sum;
    }
  }
  return dist;
}

struct W1_centr_centr_parallel : public Worker {
  NumericMatrix C;
  BMAcc_RW<double> output;
  size_t K, n;
  W1_centr_centr_parallel(NumericMatrix C, BMAcc_RW<double> output)
    : C(C),output(output) {
    K = C.ncol();
    n = C.nrow();
  }
  void operator()(size_t begin, size_t end) {
    for(size_t i = begin; i< end; i++) {
      double last = C(n-1, i);
     for(size_t j =0; j< i; j++){
       double dist = 0;

       double last_current =  C(n-1, j);

       for( size_t k =0; k < n; k++){
         dist += std::abs(C(k, j)/last_current - C(k, i)/last);
       }
        output(i, j) = dist;
       output(j,i) = dist;
     }

    }
  }
};
//' Parallel computation of the Wasserstein distances between cluster centroids
//'  
//' @param C_cumsum is a Filebacked data cummulative sum matrix of N cluster centroids
//' @param Dist is a Filebacked inter distance matrix
//' @export
//' @return N x N matrix of the Wasserstein distances
//' @keywords internal
// [[Rcpp::export]]
void W1_centr_centr(NumericMatrix C_cumsum, Environment Dist){
  XPtr<FBM_RW> xp_Dist= Dist["address_rw"];
  BMAcc_RW<double> mat_Dist(xp_Dist);
  size_t K = C_cumsum.ncol();
  W1_centr_centr_parallel DistInter(C_cumsum, mat_Dist);
  parallelFor(0, K, DistInter);

}
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

struct W1_cl_centr_parallel_BIG : public Worker {
  SubBMAcc_RW<double> M_input;
  NumericVector output;
  BMAcc_RW<double> C_cumsum;
  size_t n;
  double last;
  size_t nbr_cluster;
  W1_cl_centr_parallel_BIG(SubBMAcc_RW<double> M_input, NumericVector output, BMAcc_RW<double> C_cumsum, 
                           size_t nbr_cluster)
    : M_input(M_input),output(output), C_cumsum(C_cumsum), nbr_cluster(nbr_cluster) {
    n = M_input.nrow();

  }
  void operator()(size_t begin, size_t end) {
    for(size_t i = begin; i< end; i++) {

      double dist = 0;
      double last_current = M_input(n-1, i);
      last = C_cumsum(n-1, nbr_cluster);
      for( size_t k =0; k < n; k++){
        dist += std::abs(M_input(k, i)/last_current - C_cumsum(k,nbr_cluster)/last);
      }
      output(i) = dist;
    }
  }
};
//'@title W1_cl_centr_BIG
//'@description Parallel computation for big clusters of
//'  the Wasserstein distances between cluster centroid and signals in the cluster
//' @param X_cumsum is a Filebacked data cummulative sum matrix n x N
//' @param rowInd is a vector of the row indeces
//' @param colInd is a vector of the column indeces
//' @param C_cumsum is a vector of the cummulative sum of the cluster centroid
//' @param nbr_cluster is a number of clusters
//' @export
//' @return A vector of the Wasserstein distances
//' @keywords internal
// [[Rcpp::export]]
NumericVector W1_cl_centr_BIG(Environment X_cumsum,const IntegerVector& rowInd, const IntegerVector& colInd,
                          Environment C_cumsum, size_t nbr_cluster){

  XPtr<FBM_RW> xp_X= X_cumsum["address_rw"];
  SubBMAcc_RW<double> mat_X(xp_X, rowInd -1, colInd-1);

  XPtr<FBM_RW> xp_C_cumsum= C_cumsum["address_rw"];
  BMAcc_RW<double> mat_C_cumsum(xp_C_cumsum);
  size_t size_cl = colInd.size();
  NumericVector output(size_cl);

  W1_cl_centr_parallel_BIG Dist(mat_X, output, mat_C_cumsum, nbr_cluster);
  parallelFor(0, size_cl, Dist);
  return output;
}


struct W1_centr_centr_parallel_BIG : public Worker {
  BMAcc_RW<double> C;
  BMAcc_RW<double> output;
  size_t K, n;
  W1_centr_centr_parallel_BIG(BMAcc_RW<double> C, BMAcc_RW<double> output)
    : C(C),output(output) {
    K = C.ncol();
    n = C.nrow();
  }
  void operator()(size_t begin, size_t end) {
    for(size_t i = begin; i< end; i++) {
      double last = C(n-1, i);
      for(size_t j =0; j< i; j++){
        double dist = 0;

        double last_current =  C(n-1, j);

        for( size_t k =0; k < n; k++){
          dist += std::abs(C(k, j)/last_current - C(k, i)/last);
        }
         output(i, j) = dist;
         output(j,i) = dist;
      }

    }
  }
};

//' Parallel computation for big clusters of the Wasserstein distances between cluster centroids
//'   
//' @param C_cumsum is a Filebacked data cummulative sum matrix n x N
//' @param Dist is a Filebacked matrix of the Wasserstein distances between centroids
//' @export
//' @return NULL
//' @keywords internal
// [[Rcpp::export]]
void W1_centr_centr_BIG(Environment C_cumsum, Environment Dist){
  XPtr<FBM_RW> xp_Dist= Dist["address_rw"];
  BMAcc_RW<double> mat_Dist(xp_Dist);
  XPtr<FBM_RW> xp_C_cumsum= C_cumsum["address_rw"];
  BMAcc_RW<double> mat_C_cumsum(xp_C_cumsum);
  size_t K = mat_C_cumsum.ncol();


  W1_centr_centr_parallel_BIG DistInter(mat_C_cumsum, mat_Dist);
  parallelFor(0, K, DistInter);

}
//' @title DBindex
//' @description Davies Bouldin index computation. 
//' @param X_cumsum is a Filebacked Big Matrix, which contains column matrix cumulative sums.
//' @param rowInd is a vector of the row indeces
//' @param Clusters is a list of cluster assignment
//' @param C is a matrix of the cluster centroids
//' @param cluster_size is a vector of the cluster sizes
//' @param DistInter is a Filebacked Big Matrix, which stores distances between a cluster centoid and signals in the cluster
//' @export
//' @return A value of the DB index
//' @keywords internal
// [[Rcpp::export]]
List DBindex(Environment X_cumsum, const IntegerVector& rowInd, List Clusters, NumericMatrix C,
             const IntegerVector& cluster_size, Environment DistInter) {

  size_t K = C.ncol();
  size_t n = C.nrow();

  NumericMatrix C_cumsum(n,K);
  NumericVector dist_intra(K);
  NumericMatrix R(K,K);

  C_cumsum = cumsum_Mat(C);

  for(size_t k =0; k< K; k++){
    dist_intra(k) = sum(W1_cl_centr(X_cumsum, rowInd, Clusters[k], C_cumsum(_,k)))/cluster_size[k];

  }
   W1_centr_centr(C_cumsum, DistInter);

    XPtr<FBM_RW> xp_Dist= DistInter["address_rw"];
    BMAcc_RW<double> mat_Dist(xp_Dist);

    for(size_t i=0; i< K; i++){
      for(size_t j=0; j< i;j++){
        R(i,j) = (dist_intra(i) + dist_intra(j))/mat_Dist(i,j);
        R(j,i) = R(i,j);
      }
    }

  List res;
  res["intra"] = dist_intra;
  res["R"] = R;
  return res;
}
//' @title DBindex_BIG
//' @description Davies Bouldin index computation 
//' @param X_cumsum is a Filebacked data cummulative sum matrix n x N
//' @param rowInd is a vector of the row indeces
//' @param Clusters is a list of cluster assignment
//' @param C is a matrix of the cluster centroids
//' @param C_cumsum is a matrix of the cummulative sums of cluster centroids
//' @param cluster_size is a vector of the cluster sizes
//' @param DistInter is a matrix of the inter distances
//' @param R is a matrix of the intra distances
//' @export
//' @return A value of the DB index
//' @keywords internal
// [[Rcpp::export]]
void DBindex_BIG(Environment X_cumsum, const IntegerVector& rowInd, List Clusters, Environment C, 
                 Environment C_cumsum,
             const IntegerVector& cluster_size,  Environment DistInter, Environment R) {

  Rcout<< "Cumsum C" << endl;
  cumsum_parallel(C, C_cumsum);

  size_t K = cluster_size.size();
  NumericVector dist_intra(K);


  for(size_t k =0; k< K; k++){

    dist_intra(k) = sum(W1_cl_centr_BIG(X_cumsum, rowInd, Clusters[k], C_cumsum, k))/cluster_size[k];
//Rcout<< "Cluster"<< k<< "dist"<< dist_intra(k)<< endl;
  }

   W1_centr_centr_BIG(C_cumsum, DistInter);

    XPtr<FBM_RW> xp_Dist= DistInter["address_rw"];
    BMAcc_RW<double> mat_Dist(xp_Dist);

    XPtr<FBM_RW> xp_R= R["address_rw"];
    BMAcc_RW<double> mat_R(xp_R);
    Rcout<< "R matrix"<< endl;
    for(size_t i=0; i< K; i++){
      for(size_t j=0; j< i;j++){
         mat_R(i,j) = (dist_intra(i) + dist_intra(j))/mat_Dist(i,j);
         mat_R(j,i) = mat_R(i,j);
      }
    }
}

